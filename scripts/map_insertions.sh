#!/usr/bin/bash

# Script to map insertions from a sequencing run by BLASTing the RE of the 
# insertion sequence.
# USAGE:
#   map_insertions.sh -b brna_plasmid_file.gbk -s file1.gbk [file2.gbk ...] -o /path/to/wk_dir -f /path/to/fastq
# where brna_plasmid_file.gbk is the path to the plasmid with the bRNA 
# And file1, file2, ... are all other sequences present in the experiment

# Initialize variables
brna_plasmid_file=""
other_sequence_files=()
wk_dir=""
fastq_path=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b)
            brna_plasmid_file="$2"
            shift 2
            ;;
        -s)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                other_sequence_files+=("$1")
                shift
            done
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
if [[ -z "$brna_plasmid_file" || ${#other_sequence_files[@]} -eq 0 || -z "$wk_dir" || -z "$fastq_path" ]]; then
    echo "Usage: $0 -b brna_plasmid_file.gbk -s file1.gbk [file2.gbk ...] -o /path/to/wk_dir -f /path/to/fastq"
    exit 1
fi

mkdir -p "$wk_dir/blast" "$wk_dir/input"

python extract_insertion_end.py \
    -b "$brna_plasmid_file" \
    -s "${other_sequence_files[@]}" \
    -o "$wk_dir/input"

re_file="$wk_dir/input/RE.fasta"
genome_file="$wk_dir/input/sequences.fasta"

# Store all sequences as a single fasta
python genbank_to_fasta.py \
    "$genome_file" \
    "$brna_plasmid_file" \
    "${other_sequence_files[@]}"

# PROCESS THE SEQUENCING DATA
# if [[ "$fastq_path" == *.fastq.gz ]]; then
#     # Get the base name without .fastq.gz
#     sequencing_basename=$(basename "$fastq_path" .fastq.gz)
#     gunzip -c "$fastq_path" | seqtk seq -a > "$wk_dir/blast/$sequencing_basename.fasta"

# elif [[ "$fastq_path" == *.fastq ]]; then
#     # Get the base name without .fastq
#     sequencing_basename=$(basename "$fastq_path" .fastq)
#     seqtk seq -a "$fastq_path" > "$wk_dir/blast/$sequencing_basename.fasta"
# else
#     echo "Error: Unsupported file extension. Expected .fastq or .fastq.gz"
#     exit 1
# fi

cutadapt -g "file:$re_file" --revcomp --fasta --discard-untrimmed -O 20 -j 0 "$fastq_path" > "$wk_dir/RE_contain_contigs.fasta"
# #Header names are too long, make it shorter
# seqkit replace \
#   -p '^([^ ]+).*$' \
#   -r '$1' \
#   "$wk_dir/blast/$sequencing_basename.fasta" \
# > "$wk_dir/blast/$sequencing_basename.idonly.fasta"

# #delete temp file
# rm "$wk_dir/blast/$sequencing_basename.fasta"
# mv "$wk_dir/blast/$sequencing_basename.idonly.fasta" "$wk_dir/blast/$sequencing_basename.fasta"

# #make blast database
# makeblastdb -parse_seqids -in "$wk_dir/blast/$sequencing_basename.fasta" -dbtype nucl -out "$wk_dir/blast/$sequencing_basename"

# blastn \
# -db "$wk_dir/blast/$sequencing_basename" \
# -task blastn-short \
# -query "$re_file" \
# -perc_identity 90 -qcov_hsp_perc 90 \
# -outfmt "6 sseqid sstart send sstrand" \
# -word_size 4 \
# -soft_masking false \
# -dust no \
# -max_target_seqs 1000000 \
# -out "$wk_dir/blast/RE_blast.tsv"

# #Extract the downstream of RE from all RE-contain contigs
# python RE_blast_contig_fetch.py \
#   -f "$wk_dir/blast/$sequencing_basename.fasta" \
#   -t "$wk_dir/blast/RE_blast.tsv" \
#   -o "$wk_dir/RE_contain_contigs.fasta"
# #   -s 200 \

echo "Mapping RE onto reference genome"
minimap2 -t 48 -ax map-ont "$genome_file" "$wk_dir/RE_contain_contigs.fasta" | samtools view -bS | samtools sort -o "$wk_dir/aligned_contigs.bam" - && samtools index "$wk_dir/aligned_contigs.bam"

# Map the raw reads onto the genome, then compute coverage
echo "Mapping raw fastq reads onto reference genome"
minimap2 -t 48 -ax map-ont "$genome_file" "$fastq_path" | samtools view -bS | samtools sort -o "$wk_dir/raw_reads_aligned.bam" - && samtools index "$wk_dir/raw_reads_aligned.bam"
samtools depth -o "$wk_dir/raw_reads_coverage.tsv" "$wk_dir/raw_reads_aligned.bam"
