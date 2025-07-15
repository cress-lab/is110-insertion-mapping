import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import argparse
import os
from collections import Counter

def plot_insertions(counts_df, total_genome_length, target_coordinate=None, mode="mb", ax=None,):
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 6))
    count_df = counts_df.copy()
    # Plot Mb if mb specified as mode; else default to plotting bp
    if mode.lower() == "mb":
        position_plot_col = "position_mb"
        total_genome_length = total_genome_length / 1_000_000
    elif mode.lower() == "b":
        position_plot_col = "position"
    else:
        raise ValueError("Invalid mode for plotting! Options are mb for Mb scale or b for bp scale")
    # Set the count threshold for displaying labels.
    # A threshold of 0 will label every single point.
    label_threshold = 0
    # Filter the data to select points that meet the label threshold.
    labels_df = count_df[count_df['count'] >= label_threshold]
    print(f"  - Found {len(labels_df)} sites with counts >= {label_threshold} to label.")
    markerline, stemlines, baseline = ax.stem(
        count_df[position_plot_col],
        count_df['count'],
        linefmt='grey',
        markerfmt='o',
        basefmt=' '
    )
    plt.setp(markerline, 'markersize', 3, 'color', 'royalblue')

    # Add labels to stem plot
    for index, row in labels_df.iterrows():
        ax.text(
            x=row[position_plot_col] + 0.01,
            y=row['count'],
            s=f"{int(row['position'])}, {int(row['count'])}times",
            color='blue',
            fontsize=9
        )


    # Title and labels for ax
    ax.set_ylabel("Insertion Counts", fontsize=12)
    ax.set_ylim(bottom=0)

    # Optional x-axis limit
    if total_genome_length > 0:
        ax.set_xlim(left=0, right=total_genome_length)
    else:
        ax.set_xlim(left=0)
    return ax

def plot_triangle(target_coordinate, ax):
    # Plot a red triangle at the target site
    x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
    triangle_width = x_range * 0.01
    # triangle_height = ax.get_ylim()[1] * 0.02  
    triangle_height = ax.get_ylim()[1] * 0.022
    # print(triangle_width, triangle_height, x_range)

    triangle = patches.Polygon(
        [[target_coordinate - triangle_width / 2, -triangle_height],
        [target_coordinate + triangle_width / 2, -triangle_height],
        [target_coordinate, 0]],
        transform=ax.transData,
        clip_on=False,
        color='red',
        alpha=1
    )
    ax.add_patch(triangle)

def main():
    # Set a professional style for the plots.
    sns.set_style("whitegrid")

    # --- Step 2: Parameter Configuration ---
    # Specify the path to the SAM/BAM file for analysis.

    parser = argparse.ArgumentParser(
        description="Script to plot insertion map results against a reference genome"
    )
    parser.add_argument(
        "-d", "--wk_dir", required=True,
        help="Working directory path. Must contain a file called aligned_contigs.bam"
    )
    # parser.add_argument(
    #     "-o", "--out_dir", required=True,
    #     help="Directory to write the output plots to"
    # )
    parser.add_argument(
        "-t", "--title", required=False,
        help="Title for plot"
    )
    parser.add_argument(
        "-x", "--target_coordinate", required=False, type=int,
        help="Target coordinate for integrations in bp"
    )
    parser.add_argument(
        "-w", "--zoom_window_size", required=False, type=int, default=100,
        help="If target coordinate provided: window size for additional zoomed in plot around target coordinate"
    )
    parser.add_argument(
        "--on_target_threshold", required=False, type=int, default=10,
        help="If target coordinate provided: how many bp away from target site on either side to still consider as on target"
    )
    parser.add_argument(
        "--plot_coverage",
        action="store_true",
        help="Include this flag to enable coverage plotting. Requires a file called raw_reads_coverage.tsv in wk_dir"
    )

    args = parser.parse_args()
    alignment_file = os.path.join(args.wk_dir, "aligned_contigs.bam")
    assert os.path.exists(alignment_file), "Working dir must contain aligned_contigs.bam!"
    out_dir = os.path.join(args.wk_dir, "plots")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Set default title if not provided
    if not args.title:
        args.title = f"Insertion map for {alignment_file}"

    print(f"Step 2: Set alignment file to '{alignment_file}'")

    # --- Step 3: Data Processing (Extracting Insertion Sites and Genome Length) ---
    print("\nStep 3: Processing alignment file...")

    # Counter for tracking filtered reads.
    filter_counts = Counter({
        "unmapped": 0,
        "low_mapq": 0,
        "secondary": 0,
        "supplementary": 0,
    })
    insertion_sites, reference_names = [], []
    total_genome_length = 0

    # Open the SAM/BAM file once to process reads and get header info.
    with pysam.AlignmentFile(alignment_file, "r") as samfile:
        # Get total genome length from the SAM/BAM header.
        if samfile.header.lengths:
            total_genome_length = sum(samfile.header.lengths)
            print(f"Total genome length from header: {total_genome_length:,} bp")

        # Process all aligned reads in the file.
        for read in samfile.fetch(until_eof=True):
            # Check which filtering condition the read meets.
            if read.is_unmapped:
                filter_counts["unmapped"] += 1
                continue
            # MAPQ=0 often means a read maps to multiple locations.
            if read.mapping_quality == 0:
                filter_counts["low_mapq"] += 1
                continue
            # Secondary/supplementary alignments are alternative mappings for the same read.
            if read.is_secondary:
                filter_counts["secondary"] += 1
                continue
            if read.is_supplementary:
                filter_counts["supplementary"] += 1
                continue

            # For reads that pass the filters, get the 5' end position.
            pos = read.reference_end - 1 if read.is_reverse else read.reference_start
            insertion_sites.append(pos)
            reference_names.append(read.reference_name)

    # Report the filtering statistics.
    total_reads = sum(filter_counts.values()) + len(insertion_sites)
    print(f"Total reads processed: {total_reads:,}")
    for reason, count in filter_counts.items():
        print(f"  - {reason:14s}: {count:,} reads filtered out")
    print(f"  - Kept for analysis: {len(insertion_sites):,} reads")


    # --- Step 4: Data Aggregation and Formatting ---
    print("\nStep 4: Aggregating data...")

    positions_df = pd.DataFrame([insertion_sites, reference_names], index=['position', 'reference']).T

    # Convert to 1-based indexing
    positions_df["position"] += 1

    # Convert the results into a DataFrame for plotting.
    counts_df = pd.DataFrame(positions_df.groupby("reference")["position"].value_counts()).reset_index()

    # Add a column for plotting in Megabases (Mb).
    counts_df['position_mb'] = counts_df['position'] / 1_000_000
    # NOTE: The log10_count column is calculated but not used in the plot, so it's removed for clarity.
    print("  - Data aggregation complete.")

    # --- Step 5: Visualization (Labeling Peaks) ---
    print("\nStep 5: Creating genome-wide plot...")
    ##### Make the high level insertion map
    ax = plot_insertions(counts_df, total_genome_length, target_coordinate=args.target_coordinate, mode="mb")
    plot_triangle(args.target_coordinate / 1_000_000, ax)
    ax.set_title(args.title, fontsize=16)
    ax.set_xlabel("Position in Genome (Mb)", fontsize=12)
    plt.savefig(os.path.join(out_dir, "alignment_plot.png"), dpi=300)
    ##### If you want to plot coverage below insertion map, this will do that
    if args.plot_coverage:
        print("Plotting insertion map with coverage...")
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        plot_insertions(counts_df, total_genome_length, target_coordinate=args.target_coordinate, mode="mb", ax=ax1)
        plot_triangle(args.target_coordinate / 1_000_000, ax1) # Plot triangle at target site

        # Read dataframe of whole genome coverage
        cov_df = pd.read_csv((os.path.join(args.wk_dir, "raw_reads_coverage.tsv")), sep="\t", header=None)
        cov_df.columns=["chrom", "pos", "depth"]
        cov_df["pos_mb"] = cov_df["pos"] / 1_000_000
        # Plot coverage below insetion map
        ax2.fill_between(cov_df['pos_mb'], cov_df['depth'], color='royalblue')
        ax2.set_ylabel("Read Count", fontsize=12)
        ax2.set_xlabel("Position in Genome (Mb)", fontsize=12)
        ax2.set_ylim(bottom=0)
        ax2.set_xlim(cov_df["pos_mb"].min(), cov_df["pos_mb"].max())
        ax2.tick_params(
            direction='out',
            length=4,
            width=1,
            colors='black',
            grid_color='none',
            grid_alpha=0.5
        )
        sns.despine(ax=ax2)

        # --- Layout adjustment ---
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.1)  # Increase vertical space between plots
        plt.suptitle(args.title, fontsize=16, y=1.03)
        fig.savefig(os.path.join(out_dir, "alignment_plot_with_coverage.png"), dpi=300)
    # If target coordinate provided: save an additional zoomed in plot around 
    # that site AND save a pie chart of on target % (on target is defined as 
    # give or take args.on_target_threshold from target site [defaults to 10bp])
    
    if args.target_coordinate:
        print("Plotting zoomed in insertion map...")
        # PLOT A ZOOMED IN VERSION OF THE INSERTION MAP
        w = args.zoom_window_size / 2
        xmin, xmax = args.target_coordinate - w, args.target_coordinate + w

        ax = plot_insertions(
            counts_df[(counts_df["position"] >= xmin) & (counts_df["position"] <= xmax)], 
            total_genome_length,
            target_coordinate=args.target_coordinate, mode="b")
        plt.xlim(xmin, xmax)
        plot_triangle(args.target_coordinate, ax)
        plt.title(args.title, fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "alignment_plot_zoomed.png"), dpi=300)

        # PLOT A PIE CHART SHOWING % OF ON TARGET INTEGRATIONS
        print("Plotting pie chart...")
        def on_target_mapper(pos):
            if (pos >= (args.target_coordinate - args.on_target_threshold)) and (pos <= args.target_coordinate + args.on_target_threshold ):
                return True
            else:
                return False

        counts_df["on_target"] = counts_df["position"].map(on_target_mapper)
        counts = counts_df.groupby('on_target')['count'].sum()
        # Plot pie chart
        plt.figure(figsize=(6, 6))
        counts.plot.pie(
            autopct='%1.1f%%',
            startangle=90,
            labels=['Off Target', 'On Target'],  # Map False to Off, True to On
            colors=['lightcoral', 'lightgreen'],
            # wedgeprops={'edgecolor': 'black'}
        )
        plt.ylabel("")
        plt.savefig(os.path.join(out_dir, "on_target_pie_chart.png"), dpi=300)

if __name__ == "__main__":
    main()