"""
Script to process results of BLASTn on insertion-containing reads.
Usage: python process_16s_blastn.py <wk_dir>
Usage: python process_16s_blastn.py <out_dir> <blastn_file> <reads_file>

Outputs the following to the same directory as the blastn file: 
- 16S_reads.fasta: Reads with BLASTn hit to 16S
- reads_no_16S_blastn_hit.fasta: Reads with no hits to 16S
- 16s_blastn_hit_contig_ids.tsv: Unique contigs in 16S file that each 16S read 
  maps to
- 16s_target_seqs.tsv: Target site positions/sequences extracted from 16S 
"""
from Bio import SeqIO
import pandas as pd
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--blast_file",
    help="Path to BLAST output tsv"
)
parser.add_argument(
    "--fasta_16s",
    help="Path to 16S fasta file"
)
parser.add_argument(
    "--tseqlen",
    default=10, type=int,
    help="Path to 16S fasta file"
)
parser.add_argument(
    "--reads_file",
    help="Path to fasta with insertion-containing reads"
)
args = parser.parse_args()

OUT_DIR = os.path.dirname(args.blast_file)
# BLASTN_FILE = sys.argv[2] 
# READS_FILE = sys.argv[3] 
# Output files
summary_table_path = os.path.join(OUT_DIR, "16s_blastn_hit_contig_ids.tsv")
reads_16s_fasta = os.path.join(OUT_DIR, "16S_reads.fasta")
reads_non_16s_fasta = os.path.join(OUT_DIR, "reads_no_16S_blastn_hit.fasta")

df = pd.read_csv(args.blast_file, sep="\t", header=None, comment="#")
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
records = SeqIO.parse(args.reads_file, "fasta")
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
# Use the 16S fasta to extract target sites
rrna_16s_sequences = [r for r in SeqIO.parse(args.fasta_16s, "fasta")]
results = []
for read in df["qaccver"].unique():
    subdf = df[df["qaccver"] == read]
    # NOTE: If multiple best hits with same bitscore, only takes first
    best_hit = subdf.loc[subdf["bitscore"].idxmax()]
    best_hit_16s_name = best_hit["saccver"]
    end_position = best_hit["send"]
    for record in rrna_16s_sequences:
        if record.id == best_hit_16s_name:
            tseq = str(record.seq[end_position-args.tseqlen:end_position+args.tseqlen])
            results.append([best_hit_16s_name, end_position, "+", tseq])

rdf = pd.DataFrame(results, columns=["reference", "position", "strand", "target_seq"])
# Count unique sites
rdf = rdf.groupby(rdf.columns.tolist()).size().reset_index(name='count')
rdf.drop_duplicates(inplace=True)
# Reorder columns
rdf = rdf[["reference", "position", "strand", "count", "target_seq"]]
print(f"SAVING TO: {os.path.join(OUT_DIR, '16s_target_seqs.tsv')}")
rdf.to_csv(
    os.path.join(OUT_DIR, "16s_target_seqs.tsv"), sep="\t", index=False)


