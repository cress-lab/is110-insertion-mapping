import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from collections import Counter

# Set a professional style for the plots.
sns.set_style("whitegrid")

# --- Step 2: Parameter Configuration ---
# Specify the path to the SAM/BAM file for analysis.
# args.alignment_file = "/groups/cress/projects/efactors/ryuichi/is110/20250707/PV7SNC_3_dIS110-conj1-pool/aligned_contigs.sam"

parser = argparse.ArgumentParser(
    description="Script to plot insertion map results against a reference genome"
)
parser.add_argument(
    "-f", "--alignment_file", required=True,
    help="Path to the read alignment file"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Path to write the output plot"
)
parser.add_argument(
    "-t", "--title", required=False,
    help="Title for plot"
)

args = parser.parse_args()

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
plt.savefig(args.output, dpi=300)
plt.show()