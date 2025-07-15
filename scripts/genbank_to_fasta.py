#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def genbank_to_fasta(output_file, input_files):
    with open(output_file, 'w') as fasta_out:
        for gbk_file in input_files:
            base_name = os.path.splitext(os.path.basename(gbk_file))[0]
            for i, record in enumerate(SeqIO.parse(gbk_file, "genbank")):
                # Create a new record with the base filename as the header
                new_id = f"{base_name}_{i+1}" if i > 0 else base_name
                new_record = SeqRecord(
                    seq=record.seq,
                    id=new_id,
                    description=""
                )
                SeqIO.write(new_record, fasta_out, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py output.fasta input1.gbk [input2.gbk ...]")
        sys.exit(1)

    output_filename = sys.argv[1]
    input_filenames = sys.argv[2:]

    genbank_to_fasta(output_filename, input_filenames)
