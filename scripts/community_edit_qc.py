"""
QC stuff for community editing pipeline. 
Inputs:
    - wk_dir: location of community editing outputs
Outputs:
    - wk_dir/qc:
        - reads_no_le.fasta: reads without the LE
        - summary.txt
"""
import os
import sys
from Bio import SeqIO

def count_reads_at_each_step():
    pass

def fetch_reads_no_le():
    pass


wk_dir = sys.argv[1]
# raw_reads_file = sys.argv[2]
# # NOTE: might not handle fastq.gz...
# raw_reads = [i for i in SeqIO.parse(raw_reads_file, "fastq")]
# Make QC dir
os.makedirs(os.path.join(wk_dir, "qc"), exist_ok=True)
# Failures
reads_fail_size_filters = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "reads_fail_size_filters.fastq"), "fastq")
]

borderline = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "borderline.fasta"), "fasta")
]

chimeras = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "chimeras.fasta"), "fasta")
]

# Successes
reads_size_filtered = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "reads_size_filtered.fastq"), "fastq")
]

le_containing_amplicons_any = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "LE_containing_amplicons_any.fasta"), "fasta")
]
le_containing_amplicons = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "LE_containing_amplicons.fasta"), "fasta")
]

nonchimeras = [
    i for i in SeqIO.parse(
        os.path.join(wk_dir, "nonchimeras.fasta"), "fasta")
]
# fetch reads without LE
le_readids = [record.id for record in le_containing_amplicons]
reads_no_le = [
    record for record in reads_size_filtered if record.id not in le_readids
    ]
# fetch FULL reads that contain LE
reads_le_full = [
    record for record in reads_size_filtered if record.id in le_readids
    ]
with open(os.path.join(wk_dir, "qc", "reads_without_le.fasta"), "w") as f:
    SeqIO.write(reads_no_le, f, "fasta")
with open(os.path.join(wk_dir, "qc", "reads_with_le_untrimmed.fasta"), "w") as f:
    SeqIO.write(reads_le_full, f, "fasta")

# Now, process all the files and save a log
total_reads = len(reads_size_filtered) + len(reads_fail_size_filters)
size_filter = f"{len(reads_size_filtered)} / {total_reads} reads pass size filters"
le_filter_any = f"{len(le_containing_amplicons_any)} / {len(reads_size_filtered)} reads contain LE"
le_filter = f"{len(le_containing_amplicons)} / {len(reads_size_filtered)} reads contain only 1 LE"
# Chimera stats
chimera_lengths = [len(record.seq) for record in chimeras]
nonchimera_lengths = [len(record.seq) for record in nonchimeras]
borderline_lengths = [len(record.seq) for record in borderline]
total_checked = len(chimeras) + len(nonchimeras) + len(borderline)
chimera_filter = (
    f"{len(nonchimeras)} nonchimeras ({len(nonchimeras)/total_checked:.1%}), "
    f"{len(chimeras)} chimeras ({len(chimeras)/total_checked:.1%}), "
    f"{len(borderline)} borderline ({len(borderline)/total_checked:.1%})"
)

summary_lines = [
    size_filter,
    le_filter_any,
    le_filter,
    chimera_filter
]

with open(os.path.join(wk_dir, "qc", "summary.txt"), "w") as f:
    for line in summary_lines:
        f.write(line + "\n")