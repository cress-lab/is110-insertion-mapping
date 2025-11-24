import pysam
import pandas as pd
from collections import Counter
import sys


def check_no_overlap(*args):
    """
    Ensures that there is no overlap between any of the provided lists.
    """
    all_items = set()
    for lst in args:
        for item in lst:
            if item in all_items:
                return False
            all_items.add(item)
    return True


def parse_bam_mapping_locations(alignment_file):
    with pysam.AlignmentFile(alignment_file, "r") as samfile:
        # filter_counts = Counter({
        #     "unmapped": 0,
        #     "low_mapq": 0,
        #     "secondary": 0,
        #     "supplementary": 0,
        # })
        unmapped_reads, ambiguous_mapped_reads, mapped_reads = set(), set(), set()
        mapping_locations = dict() # Store all locations each read maps to

        if samfile.header.lengths:
            total_genome_lengths = dict(
                zip(samfile.header.references, samfile.header.lengths))
            print(f"Total genome lengths from header:")
            for genome, length in total_genome_lengths.items():
                print(f"\t{genome}: {length:,} bp")
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped:
                unmapped_reads.add(read.query_name)
                continue
            elif read.mapping_quality == 0:
                # Get read alignment quality
                # filter_counts["low_mapq"] += 1
                if (not read.is_secondary) and (not read.is_supplementary):
                    ambiguous_mapped_reads.add(read.query_name)
            else:
                if read.is_secondary or read.is_supplementary:
                    continue
                mapped_reads.add(read.query_name)
            # Save where the read is mapped to
            mapped_to = read.reference_name
            if read.query_name not in mapping_locations:
                mapping_locations[read.query_name] = [mapped_to]
            else:
                if mapped_to not in mapping_locations[read.query_name]:
                    mapping_locations[read.query_name].append(mapped_to)

        # convert to lists
        unmapped_reads = list(unmapped_reads)
        ambiguous_mapped_reads = list(ambiguous_mapped_reads)
        mapped_reads = list(mapped_reads)
        assert check_no_overlap(
            unmapped_reads, ambiguous_mapped_reads, mapped_reads
            ), "Overlap found between read categories!"
        # Get a dataframe of mapping locations
        for k, v in mapping_locations.items():
            mapping_locations[k] = ";".join(v)
        map_location_df = pd.DataFrame.from_dict(
            mapping_locations, orient='index', columns=['mapped_to']).reset_index().rename(columns={'index': 'read_name'}) 
        return unmapped_reads, ambiguous_mapped_reads, mapped_reads, map_location_df


if __name__ == "__main__":
    WK_DIR = sys.argv[1]
    bam_16s = parse_bam_mapping_locations(f"{WK_DIR}/mapping_16s/16S_map.bam")
    non_16s, ambiguous_reads_16s, mapped_reads_16s, mapping_locations_16s = bam_16s
    reads_16s = ambiguous_reads_16s + mapped_reads_16s
    bam_genomes = parse_bam_mapping_locations(f"{WK_DIR}/mapping_16s/non_16S_mapped_to_genomes.bam")
    unmapped, multimapped_to_genome, mapped_to_genome, _ = bam_genomes
    counts_dict = {
        "Category": [
            "16S reads",
            "Non-16S",
            "Unmapped to genome"
        ],
        "Count": [
            len(reads_16s),
            len(mapped_to_genome + multimapped_to_genome),
            len(unmapped)
        ]
    }

    counts_df = pd.DataFrame(counts_dict)
    counts_df.to_csv(f"{WK_DIR}/mapping_16s/16S_mapping_summary.tsv", index=False, sep="\t")
    mapping_locations_16s.to_csv(f"{WK_DIR}/mapping_16s/16S_mapping_locations.tsv", index=False, sep="\t")
