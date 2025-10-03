#!/usr/bin/bash

# Script to determine taxonomy for community editing amplicons. 
# USAGE:
#   ./community_editing.sh -b brna_plasmid_file.gbk -o /path/to/wk_dir -f /path/to/fastq -n max_length
# where brna_plasmid_file.gbk is the path to the plasmid with the bRNA 
# NOTE: Assumes there exists a file $wk_dir/LE.fasta

# Initialize variables
# brna_plasmid_file=""
# other_sequence_files=()
wk_dir=""
fastq_path=""
max_length=575

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        # -b)
        #     brna_plasmid_file="$2"
        #     shift 2
        #     ;;
        -n)
            max_length="$2"
            shift 2
            ;;
        -o)
            wk_dir="$2"
            shift 2
            ;;
        -f)
            fastq_path="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required inputs
# if [[ -z "$brna_plasmid_file" || -z "$wk_dir" || -z "$fastq_path" ]]; then
if [[ -z "$wk_dir" || -z "$fastq_path" ]]; then
    echo "Usage: $0 -b brna_plasmid_file.gbk -o /path/to/wk_dir -f /path/to/fastq"
    exit 1
fi

mkdir -p "$wk_dir/input"

# NOTE: 
# conda activate brna
# python extract_insertion_end.py \
#     -b "$brna_plasmid_file" \
#     -o "$wk_dir/input" \
#     -m "L" \
#     -n 50

python preprocess_community_editing.py \
    -i "$fastq_path" \
    -n "$max_length" \
    -o "$wk_dir"

# conda activate qiime2-amplicon-2024.10
end_file="$wk_dir/input/LE.fasta"
# PROCESS THE SEQUENCING DATA
cutadapt -a "file:$end_file" --revcomp  --fasta --discard-untrimmed  -O 10 -j 0 "$wk_dir/reads.fastq" > "$wk_dir/LE_containing_amplicons.fasta"
echo "Running analysis on sample $wk_dir..."

vsearch --uchime_ref "$wk_dir/LE_containing_amplicons.fasta" \
    --db /shared/db/silva/latest/fasta/SILVA_138_SSURef_NR99_tax_silva_trunc.replace-whitespace.fasta \
    --nonchimeras "$wk_dir/nonchimeras.fasta" \
    --chimeras "$wk_dir/chimeras.fasta" \
    --borderline "$wk_dir/borderline.fasta"

qiime tools import \
    --input-path "$wk_dir/nonchimeras.fasta" \
    --output-path "$wk_dir/nonchimeras.qza" \
    --type 'FeatureData[Sequence]' \
    --input-format DNAFASTAFormat

qiime feature-classifier classify-sklearn \
    --i-classifier /groups/cress/projects/efactors/ryuichi/is110/20250808/qiime/silva-138-99-nb-classifier.qza \
    --i-reads "$wk_dir/nonchimeras.qza" \
    --o-classification "$wk_dir/nonchimeras.taxonomy.qza"

qiime tools export \
    --input-path "$wk_dir/nonchimeras.taxonomy.qza" \
    --output-path "$wk_dir/taxonomy"

python process_taxonomy.py "$wk_dir"
