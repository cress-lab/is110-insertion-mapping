import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import argparse
import os
from collections import Counter

# Set a professional style for the plots.
sns.set_style("whitegrid")

# --- Step 2: Parameter Configuration ---
# Specify the path to the SAM/BAM file for analysis.

parser = argparse.ArgumentParser(
    description="Script to plot insertion map results against a reference genome"
)
parser.add_argument(
    "-f", "--alignment_file", required=True,
    help="Path to the read alignment file"
)
parser.add_argument(
    "-o", "--out_dir", required=True,
    help="Directory to write the output plots to"
)
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

args = parser.parse_args()

def add_triangle(x, size_frac=0.01):
    """
    Add an upward-pointing triangle below the x-axis at data coordinate x,
    with width = height, both relative to axes size.

    Parameters:
    - x: data coordinate on x-axis where triangle is centered.
    - size_frac: fraction of axes width (and height) for the triangle size.
    """
    ax = plt.gca()
    fig = plt.gcf()

    # Get data limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    axes_width_data = xlim[1] - xlim[0]
    axes_height_data = ylim[1] - ylim[0]

    # Get axes bounding box in pixels
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes_width_inch = bbox.width
    axes_height_inch = bbox.height

    # Calculate scale factors: data units per inch
    x_scale = axes_width_data / axes_width_inch
    y_scale = axes_height_data / axes_height_inch

    # size_frac of axes width in data units (x direction)
    width = axes_width_data * size_frac

    # To get equal width and height visually, height in data units is width * (y_scale / x_scale)
    height = width * (y_scale / x_scale)

    # Define transform: x in data coords, y in axis coords (0 is bottom of axes)
    transform = ax.get_xaxis_transform()

    # Convert height to axis coords:
    # Since y in transform is axis coords (0 to 1), convert height from data to axis coords:
    height_axis = height / axes_height_data  # height in axis fraction

    triangle = patches.Polygon(
        [[x - width / 2, -height_axis],
         [x + width / 2, -height_axis],
         [x, 0]],
        transform=transform,
        clip_on=False,
        color='red',
        alpha=1
    )
    ax.add_patch(triangle)
    plt.subplots_adjust(bottom=0.2)

    return triangle

# Set default title if not provided
if not args.title:
    args.title = f"Insertion map for {args.alignment_file}"

print(f"Step 2: Set alignment file to '{args.alignment_file}'")

# --- Step 3: Data Processing (Extracting Insertion Sites and Genome Length) ---
print("\nStep 3: Processing alignment file...")

# Counter for tracking filtered reads.
filter_counts = Counter({
    "unmapped": 0,
    "low_mapq": 0,
    "secondary": 0,
    "supplementary": 0,
})
insertion_sites = []
total_genome_length = 0

# Open the SAM/BAM file once to process reads and get header info.
with pysam.AlignmentFile(args.alignment_file, "r") as samfile:
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

# Report the filtering statistics.
total_reads = sum(filter_counts.values()) + len(insertion_sites)
print(f"Total reads processed: {total_reads:,}")
for reason, count in filter_counts.items():
    print(f"  - {reason:14s}: {count:,} reads filtered out")
print(f"  - Kept for analysis: {len(insertion_sites):,} reads")


# --- Step 4: Data Aggregation and Formatting ---
print("\nStep 4: Aggregating data...")

# Convert the list of insertion sites into a Pandas Series.
sites_series = pd.Series(insertion_sites, name="position")

#Convert to 1-based indexing.
sites_series += 1

# Count the occurrences of each position and sort by position.
position_counts = sites_series.value_counts().sort_index()

# (Optional) Print the counts for each insertion site.
# for pos, count in position_counts.items():
#     print(f"  - Insertion site at {pos}: {count} times")

# Convert the results into a DataFrame for plotting.
counts_df = position_counts.reset_index()
counts_df.columns = ['position', 'count']

# Add a column for plotting in Megabases (Mb).
counts_df['position_mb'] = counts_df['position'] / 1_000_000
# NOTE: The log10_count column is calculated but not used in the plot, so it's removed for clarity.
print("  - Data aggregation complete.")


# --- Step 5: Visualization (Labeling Peaks) ---
print("\nStep 5: Creating genome-wide plot...")

# Set the count threshold for displaying labels.
# A threshold of 0 will label every single point.
label_threshold = 0

# Filter the data to select points that meet the label threshold.
labels_df = counts_df[counts_df['count'] >= label_threshold]
print(f"  - Found {len(labels_df)} sites with counts >= {label_threshold} to label.")

# --- Plot Creation ---
plt.figure(figsize=(16, 6))

# Create the stem plot using raw counts (linear scale).
markerline, stemlines, baseline = plt.stem(
    counts_df['position_mb'],
    counts_df['count'],
    linefmt='grey', markerfmt='o', basefmt=' '
)
plt.setp(markerline, 'markersize', 3, 'color', 'royalblue')

# --- Labeling Process ---
# Loop through the filtered DataFrame to add labels.
for index, row in labels_df.iterrows():
    # Use plt.text() to add text to the plot.
    plt.text(
        x=row['position_mb'] + 0.01,  # x-coordinate: slightly offset to the right.
        y=row['count'],               # y-coordinate: same height as the marker.
        s=f"{int(row['position'])}, {int(row['count'])}times", # s: the string to display.
        color='blue',                 # color: text color.
        fontsize=9                    # fontsize: text size.
    )

# Add a red triangle below the target site if target_coordinate argument provided
if args.target_coordinate:
    # args.target_coordinate = args.target_coordinate / 1_000_000
    # tri = add_triangle(args.target_coordinate / 1_000_000)
    target_coordinate = args.target_coordinate / 1_000_000
    # Define x-position for the triangle
    triangle_height = 0.02  # in axis coordinates (relative to axes height)
    triangle_width = triangle_height * 1.8

    # Get current axes
    ax = plt.gca()

    # Define triangle in axis coordinates (y=0 is bottom of axes)
    transform = ax.get_xaxis_transform()  # Transforms x in data, y in axes coordinates

    # Create triangle below the x-axis
    triangle = patches.Polygon(
        [[target_coordinate - triangle_width / 2, -triangle_height],
        [target_coordinate + triangle_width / 2, -triangle_height],
        [target_coordinate, 0]],  # pointing up toward axis
        transform=transform,
        clip_on=False,  # Draw outside plot area
        color='red', alpha=1
    )

    # Add triangle
    ax.add_patch(triangle)

    # Optional: make room below axis
    plt.subplots_adjust(bottom=0.2)


# Set labels and title.
plt.title(args.title, fontsize=16)
plt.xlabel("Position in Genome (Mb)", fontsize=12)
# The y-axis label is corrected to reflect the linear scale of the data being plotted.
plt.ylabel("Insertion Counts", fontsize=12)
plt.ylim(bottom=0)

if total_genome_length > 0:
    total_genome_length_mb = total_genome_length / 1_000_000
    plt.xlim(left=0, right=total_genome_length_mb)
else:
    plt.xlim(left=0)

plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "alignment_plot.png"), dpi=300)
# # If target coordinate provided: save an additional zoomed in plot around that 
# # site AND save a pie chart of on target % (on target is defined as give or take 
# # 1bp from target site)
if args.target_coordinate:
    # PLOT A ZOOMED IN VERSION OF THE INSERTION MAP
    w = args.zoom_window_size / 2
    xmin, xmax = args.target_coordinate - w, args.target_coordinate + w

    plt.figure(figsize=(16, 6))
    plt.xlim(xmin, xmax)

    # Create the stem plot using raw counts (linear scale).
    markerline, stemlines, baseline = plt.stem(
        counts_df['position'],
        counts_df['count'],
        linefmt='grey', markerfmt='o', basefmt=' '
    )
    plt.setp(markerline, 'markersize', 3, 'color', 'royalblue')

    # --- Labeling Process ---
    # Loop through the filtered DataFrame to add labels.
    for index, row in labels_df.iterrows():
        if (row["position"] >= xmin) & (row["position"] <= xmax):
            # Use plt.text() to add text to the plot.
            plt.text(
                x=row['position'] + 0.01,  # x-coordinate: slightly offset to the right.
                y=row['count'],               # y-coordinate: same height as the marker.
                s=f"{int(row['position'])}, {int(row['count'])}times", # s: the string to display.
                color='blue',                 # color: text color.
                fontsize=9                    # fontsize: text size.
            )
    # Set labels and title.
    plt.title(args.title, fontsize=16)
    plt.xlabel("Position in Genome (b)", fontsize=12)
    # The y-axis label is corrected to reflect the linear scale of the data being plotted.
    plt.ylabel("Insertion Counts", fontsize=12)
    plt.ylim(bottom=0)
    tri = add_triangle(args.target_coordinate)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "alignment_plot_zoomed.png"), dpi=300)

    # PLOT A PIE CHART SHOWING % OF ON TARGET INTEGRATIONS
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
    plt.savefig(os.path.join(args.out_dir, "on_target_pie_chart.png"), dpi=300)

# plt.show()
