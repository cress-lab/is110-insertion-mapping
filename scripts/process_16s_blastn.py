"""
Script to process results of BLASTn on insertion-containing reads.
Usage: python process_16s_blastn.py <wk_dir>
Usage: python process_16s_blastn.py <out_dir> <blastn_file> <reads_file>

Outputs the following to out_dir: 
- 16S_reads.fasta: Reads with BLASTn hit to 16S
- reads_no_16S_blastn_hit.fasta: Reads with no hits to 16S
- 16s_blastn_hit_contig_ids.tsv: Unique contigs in 16S file that each 16S read 
  maps to
"""
from Bio import SeqIO
import pandas as pd
import sys
import os

OUT_DIR = sys.argv[1] 
BLASTN_FILE = sys.argv[2] 
READS_FILE = sys.argv[3] 
# Output files
summary_table_path = os.path.join(OUT_DIR, "16s_blastn_hit_contig_ids.tsv")
reads_16s_fasta = os.path.join(OUT_DIR, "16S_reads.fasta")
reads_non_16s_fasta = os.path.join(OUT_DIR, "reads_no_16S_blastn_hit.fasta")

df = pd.read_csv(BLASTN_FILE, sep="\t", header=None, comment="#")
df.columns = [
    "qaccver",
    "saccver",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
]
# Save 16S hits for each read
hits_per_query = (
    df.groupby("qaccver")["saccver"]
    .apply(lambda x: ";".join(x))
    .reset_index()
    )
hits_per_query.to_csv(summary_table_path, sep="\t", index=False)
print(f"Saved mapping locations to {summary_table_path}")
# Save reads with and without 16S
reads_16s_ids = df["qaccver"].unique().tolist()
records = SeqIO.parse(READS_FILE, "fasta")
reads_16S, reads_other = [], []
for record in records:
    if record.id in reads_16s_ids:
        reads_16S.append(record)
    else:
        reads_other.append(record)
print("Number of 16S records collected:", len(reads_16S))
print("Number of non-16S records collected:", len(reads_other))
with open(reads_16s_fasta, "w") as f1:
    SeqIO.write(reads_16S, f1, "fasta")
with open(reads_non_16s_fasta, "w") as f2:
    SeqIO.write(reads_other, f2, "fasta")
print(f"Saved reads to {reads_16s_fasta} and {reads_non_16s_fasta}")
