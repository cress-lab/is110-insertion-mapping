from Bio import SeqIO
from os.path import join
import argparse

def main():
    parser = argparse.ArgumentParser(
        description = ("""
        Script to filter amplicons by length. 
        """)
    )
    parser.add_argument("-i", "--input", help="Reads file (fastq format)",)
    parser.add_argument("-a", "--max_size", 
                        type=int, default=600,
                        help="Maximum size of reads to keep",)
    parser.add_argument("-b", "--min_size", 
                        type=int, default=400,
                        help="Minimum size of reads to keep",)
    parser.add_argument("-o", "--output", 
                        required=True, help="Output directory",)
    args = parser.parse_args()

    records = SeqIO.parse(args.input, "fastq")
    records_pass_filter, others = [], []
    for record in records:
        if (len(record.seq) <= args.max_size) & (len(record.seq) >= args.min_size):
            records_pass_filter.append(record)
        else:
            others.append(record)
    
    with open(join(args.output, "reads_size_filtered.fastq"), "w") as f:
        SeqIO.write(records_pass_filter, f, "fastq")
    with open(join(args.output, "reads_fail_size_filters.fastq"), "w") as f:
        SeqIO.write(others, f, "fastq")

if __name__ == "__main__":
    main()
