#!/usr/bin/bash

# Script to map insertions from a sequencing run by BLASTing the RE of the 
# insertion sequence.
# USAGE:
#   map_insertions.sh -b brna_plasmid_file.gbk -s file1.gbk [file2.gbk ...] -o /path/to/wk_dir -f /path/to/fastq
# where brna_plasmid_file.gbk is the path to the plasmid with the bRNA 
# And file1, file2, ... are all other sequences present in the experiment

# Initialize variables
brna_plasmid_file=""
# other_sequence_files=()
wk_dir=""
fastq_path=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b)
            brna_plasmid_file="$2"
            shift 2
            ;;
        # -s)
        #     shift
        #     while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
        #         other_sequence_files+=("$1")
        #         shift
        #     done
        #     ;;
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
if [[ -z "$brna_plasmid_file" || -z "$wk_dir" || -z "$fastq_path" ]]; then
    echo "Usage: $0 -b brna_plasmid_file.gbk -o /path/to/wk_dir -f /path/to/fastq"
    exit 1
fi

mkdir -p "$wk_dir/blast" "$wk_dir/input"

# NOTE: 
python extract_insertion_end.py \
    -b "$brna_plasmid_file" \
    -s "${other_sequence_files[@]}" \
    -o "$wk_dir/input" \
    -m "L" \
    -n 50

end_file="$wk_dir/input/LE.fasta"
# genome_file="$wk_dir/input/sequences.fasta"

# # Store all sequences as a single fasta
# python genbank_to_fasta.py \
#     "$genome_file" \
#     "$brna_plasmid_file" \
#     "${other_sequence_files[@]}"

# PROCESS THE SEQUENCING DATA
cutadapt -a "file:$end_file" --revcomp -o "$wk_dir/LE_containing_amplicons.fasta" --fasta --discard-untrimmed -O 10 -j 0 "$fastq_path"

# conda activate qiime2-amplicon-2024.10

wk_dir="$1"
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

# echo "Mapping RE onto reference genome"
# minimap2 -t 48 -ax map-ont "$genome_file" "$wk_dir/LE_containing_amplicons.fasta" | samtools view -bS | samtools sort -o "$wk_dir/aligned_contigs.bam" - && samtools index "$wk_dir/aligned_contigs.bam"

# # Map the raw reads onto the genome, then compute coverage
# echo "Mapping raw fastq reads onto reference genome"
# minimap2 -t 48 -ax map-ont "$genome_file" "$fastq_path" | samtools view -bS | samtools sort -o "$wk_dir/raw_reads_aligned.bam" - && samtools index "$wk_dir/raw_reads_aligned.bam"
# samtools depth -o "$wk_dir/raw_reads_coverage.tsv" "$wk_dir/raw_reads_aligned.bam"
