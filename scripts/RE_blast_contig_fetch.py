#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq

def main():
    """
    Main function to parse arguments and orchestrate sequence extraction.
    """
    parser = argparse.ArgumentParser(
        # Updated the description to reflect the new functionality
        description="From each BLAST hit, extract the downstream sequence all the way to the end of the contig."
    )
    parser.add_argument(
        "-f", "--fasta", required=True,
        help="Path to the indexed FASTA file (e.g., your genome or assembly)"
    )
    parser.add_argument(
        "-t", "--hits", required=True,
        help="Path to the BLAST hits TSV (sseqid, sstart, send, sstrand)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to write the output FASTA (containing the downstream regions)"
    )
    args = parser.parse_args()

    # Ensure the output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Use SeqIO.index for efficient, random access to FASTA sequences
    print(f"Indexing FASTA file: {args.fasta}...")
    try:
        seq_dict = SeqIO.index(args.fasta, "fasta")
    except Exception as e:
        print(f"Error: Could not index FASTA file. Please ensure it is a valid FASTA. Details: {e}")
        return

    print(f"Processing hits from {args.hits} and extracting downstream regions...")
    with open(args.hits) as hits_file, open(args.output, "w") as out_f:
        for line_num, line in enumerate(hits_file, 1):
            if not line.strip():
                continue
            
            try:
                id, sstart, send, strand = line.strip().split("\t")
                sstart, send = int(sstart), int(send)
            except ValueError:
                print(f"Warning: Skipping malformed line {line_num} in hits file: {line.strip()}")
                continue

            # Check if the contig ID from the hits file exists in the FASTA file
            if id not in seq_dict:
                print(f"Warning: Contig ID '{id}' on line {line_num} not found in FASTA file. Skipping.")
                continue

            # Get the full contig length
            contig_len = len(seq_dict[id])
            subseq = ""

            # --- LOGIC CHANGE START ---
            # Calculate coordinates to get the sequence to the end of the contig.

            if strand.lower() == "plus":
                # On the plus strand, "downstream" is from the hit's end to the contig's end.
                start_pos = send + 1
                end_pos = contig_len
                
                # Extract sequence only if the hit is not at the very end of the contig
                if start_pos <= contig_len:
                    subseq = seq_dict[id].seq[start_pos - 1 : end_pos]

            else:  # minus strand
                # On the minus strand, "downstream" (towards 3' end) is from the hit's start 
                # towards the beginning of the contig's coordinate system.
                start_pos = 1
                end_pos = sstart - 1

                # Extract sequence only if the hit is not at the very beginning of the contig
                if end_pos >= 1:
                    raw_seq = seq_dict[id].seq[start_pos - 1 : end_pos]
                    subseq = raw_seq.reverse_complement()
            
            # --- LOGIC CHANGE END ---

            # Write the sequence to the output file, but only if it's not empty.
            if subseq:
                # Create a more descriptive header for each unique extraction
                header = f">{id}|hit:{sstart}-{send}|strand:{strand}|downstream_region"
                out_f.write(f"{header}\n{str(subseq)}\n")

    print("\nProcess complete.")
    print(f"All available downstream sequences have been written to: {args.output}")


if __name__ == "__main__":
    main()