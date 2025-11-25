import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import argparse
import os
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
    description="Visualize insertion mapping results from BAM files",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
    Examples:
    # Basic usage with working directory only
    python visualize_insertions.py /path/to/workdir
    
    # Specify a single target position
    python visualize_insertions.py /path/to/workdir --target 12345
    
    # Specify multiple targets from a file
    python visualize_insertions.py /path/to/workdir --target targets.tsv
    
    # Specify target position with custom reference genome name
    python visualize_insertions.py /path/to/workdir --target 12345 --reference my_genome
        """
    )
    
    parser.add_argument(
        'workdir',
        help='Working directory containing BAM files and analysis results'
    )

    parser.add_argument(
        '--mapping_file', "-m",
        default="aligned_contigs.bam",
        help='Path to BAM file with mapped reads relative to working directory'
    )

    parser.add_argument(
        '--end', "-e",
        default="R",
        choices=["R", "L"],
        help='Which end of the insertion was used to find insertion-containing' \
        'reads. If R, positions of insertions will be the 5\' ends of reads; ' \
        'if L, will use the 3\' end.'
    )
    
    parser.add_argument(
        '--target', '-t',
        help='Target position(s) to highlight. Can be: '
             '1) A single integer position, or '
             '2) A .tsv file containing contig names and target positions (one per line)'
    )
    
    parser.add_argument(
        '--reference', '-r',
        help='Reference genome contig name within sequences file (required when ' \
        'target is a single integer)'
    )

    parser.add_argument(
        '--plot_coverage', 
        action='store_true', default=False,
        help='Option to make coverage plots. Defaults to false. '
    )
    
    args = parser.parse_args()
    # Set output directory
    args.out_dir = os.path.join(args.workdir, "plots")
    # Validate args
    if args.target and not args.target.endswith('.tsv'):
        try:
            int(args.target)
        except ValueError:
            parser.error("Target must be either an integer or a .tsv file")
        
        if not args.reference:
            parser.error("Reference genome name is required when target is a single integer")

    return args


def plot_insertions(counts_df, total_genome_length, mode="mb", ax=None,):
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
    return fig, ax

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


def parse_bam_for_plots(wk_dir, bam_file, end="R", title=""):
    alignment_file = os.path.join(wk_dir, bam_file)
    assert os.path.exists(alignment_file), f"{alignment_file} not found!"
    # Set default title if not provided
    if not title:
        title = f"Insertion map for {alignment_file}"

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
    insertion_sites, reference_names, read_names = [], [], []
    # total_genome_length = 0

    # Open the SAM/BAM file once to process reads and get header info.
    with pysam.AlignmentFile(alignment_file, "r") as samfile:
        # Get total genome length from the SAM/BAM header.
        if samfile.header.lengths:
            total_genome_lengths = dict(zip(samfile.header.references, samfile.header.lengths))
            # print(f"Total genome lengths from header:")
            # for genome, length in total_genome_lengths.items():
            #     print(f"\t{genome}: {length:,} bp")

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
            
            # For reads that pass the filters, get the 5' end position if using 
            # RE to find insertions, otherwise get 3'.
            if end == "R":
                pos = read.reference_end - 1 if read.is_reverse else read.reference_start
            elif end == "L":
                pos = read.reference_start if read.is_reverse else read.reference_end - 1
            insertion_sites.append(pos)
            reference_names.append(read.reference_name)
            read_names.append(read.query_name)

    # Report the filtering statistics.
    total_reads = sum(filter_counts.values()) + len(insertion_sites)
    print(f"Total reads processed: {total_reads:,}")
    for reason, count in filter_counts.items():
        print(f"  - {reason:14s}: {count:,} reads filtered out")
    print(f"  - Kept for analysis: {len(insertion_sites):,} reads")


    # --- Step 4: Data Aggregation and Formatting ---
    print("\nStep 4: Aggregating data...")

    positions_df = pd.DataFrame(
        [insertion_sites, reference_names, read_names], 
        index=['position', 'reference', "read"]).T
    positions_df = pd.DataFrame(
        [read_names, reference_names, insertion_sites], 
        index=['read', 'reference', 'position']).T

    # Convert to 1-based indexing
    positions_df["position"] += 1

    # Convert the results into a DataFrame for plotting.
    counts_df = pd.DataFrame(positions_df.groupby("reference")["position"].value_counts()).reset_index()

    # Add a column for plotting in Megabases (Mb).
    counts_df['position_mb'] = counts_df['position'] / 1_000_000
    # NOTE: The log10_count column is calculated but not used in the plot, so it's removed for clarity.
    print("  - Data aggregation complete.")
    return counts_df, positions_df, total_genome_lengths, title


def plot_coverage(args):
    cov_df = pd.read_csv((os.path.join(args.workdir, "raw_reads_coverage.tsv")), sep="\t", header=None)
    cov_df.columns=["chrom", "pos", "depth"]
    cov_df["pos_mb"] = cov_df["pos"] / 1_000_000
    for contig in cov_df["chrom"].unique():
        cov_df_contig = cov_df[cov_df["chrom"]==contig]
        # Plot coverage below insetion map
        fig, ax = plt.subplots(figsize=(8,2))
        ax.fill_between(cov_df_contig['pos_mb'], cov_df_contig['depth'], color='royalblue')
        ax.set_ylabel("Read Count", fontsize=12)
        ax.set_xlabel("Position in Genome (Mb)", fontsize=12)
        ax.set_ylim(bottom=0)
        ax.set_xlim(cov_df_contig["pos_mb"].min(), cov_df_contig["pos_mb"].max())
        ax.tick_params(
            direction='out',
            length=4,
            width=1,
            colors='black',
            grid_color='none',
            grid_alpha=0.5
        )
        sns.despine(ax=ax)
        ax.set_title(contig)

        # --- Layout adjustment ---
        plt.tight_layout()
        fig.savefig(os.path.join(args.out_dir, f"{contig}_coverage.png"))



def main():
    args = parse_args()
    # Loading stuff
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
    counts_df, positions_df, total_genome_lengths, title = parse_bam_for_plots(
        args.workdir, args.mapping_file, args.end)
    # Handle target info
    targets = None
    reference_genome_name = args.reference
    if args.target is not None:
        use_target_tsv = args.target.endswith(".tsv")
        if use_target_tsv:
            target_df = pd.read_csv(args.target, sep="\t")
            # reference_genome_name = target_df.columns[0]
            # targets = target_df.iloc[:, 0].astype(int).tolist()
        else:
            targets = [int(args.target)]
    print("\nStep 5: Creating genome-wide plots...")
    ##### Make the high level insertion map
    sns.set_style("whitegrid")
    for reference in counts_df["reference"].unique():
        fig, ax = plot_insertions(
            counts_df[counts_df["reference"]==reference], 
            total_genome_lengths[reference], 
            mode="b")
        if args.target is not None:
            if use_target_tsv:
                target_df_ref = target_df[target_df['contig']==reference]
                targets = target_df_ref['target_site'].astype(int).tolist()
                for t in targets:
                    plot_triangle(t, ax)
            # SINGLE INTEGER TARGET VERSION
            else:
                if reference == reference_genome_name:
                    for t in targets:
                        plot_triangle(t, ax)

        title = f"Insertion map for {reference}"
        ax.set_title(title, fontsize=16)
        ax.set_xlabel("Position in Genome (b)", fontsize=12)
        fig.savefig(os.path.join(args.out_dir, f"{reference}_insertion_map.png"))
        # plt.title()
        # plt.show()
    if args.plot_coverage:
        plot_coverage(args)
    # Save the counts and positions dataframes
    counts_df.to_csv(
        os.path.join(args.out_dir, "mapped_reads_per_position.tsv"), 
        sep="\t", 
        index=False
        )
    positions_df.to_csv(
        os.path.join(args.out_dir, "positions_per_read.tsv"), 
        sep="\t", 
        index=False
        )


if __name__ == "__main__":
    main()
