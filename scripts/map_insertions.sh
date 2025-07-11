#!/usr/bin/bash

# Script to map insertions from a sequencing run by BLASTing the RE of the 
# insertion sequence.
# USAGE:
#   ./map_insertions.sh <fastq path> <genome> <re> <working directory>
# where fastq path is the path to the sequencing results, genome is the path to 
# the genome, and re is the path to a fasta with the RE sequence

fastq_path="$1"
genome_file="$2"
re_seq="$3"
wk_dir="${4:-.}"

mkdir -p "$wk_dir/blast"

# Unzip the sequencing read file if it's compressed, if not don't do anything
if [[ "$fastq_path" == *.fastq.gz ]]; then
    # Get the base name without .fastq.gz
    sequencing_basename=$(basename "$fastq_path" .fastq.gz)

    # Unzip only if the .fastq doesn't already exist
    if [[ ! -f "${fastq_path%.gz}" ]]; then
        echo "Unzipping $fastq_path..."
        gunzip -d "$fastq_path"
    else
        echo "Unzipped file already exists: ${fastq_path%.gz}"
    fi

elif [[ "$fastq_path" == *.fastq ]]; then
    # Get the base name without .fastq
    sequencing_basename=$(basename "$fastq_path" .fastq)
else
    echo "Error: Unsupported file extension. Expected .fastq or .fastq.gz"
    exit 1
fi

# sequencing_basename=$(basename "$fastq_path" .fastq.gz)

# Convert fastq to fasta
seqtk seq -a "$wk_dir/input/$sequencing_basename.fastq" > "$wk_dir/blast/$sequencing_basename.fasta"

#Header names are too long, make it shorter
seqkit replace \
  -p '^([^ ]+).*$' \
  -r '$1' \
  "$wk_dir/blast/$sequencing_basename.fasta" \
> "$wk_dir/blast/$sequencing_basename.idonly.fasta"

#delete temp file
rm "$wk_dir/blast/$sequencing_basename.fasta"; mv "$wk_dir/blast/$sequencing_basename.idonly.fasta" "$wk_dir/blast/$sequencing_basename.fasta"

#make blast database
makeblastdb -parse_seqids -in "$wk_dir/blast/$sequencing_basename.fasta" -dbtype nucl -out "$wk_dir/blast/$sequencing_basename"

blastn \
-db "$wk_dir/blast/$sequencing_basename" \
-task blastn-short \
-query "$re_seq" \
-perc_identity 100 -qcov_hsp_perc 100 \
-outfmt "6 sseqid sstart send sstrand" \
-word_size 4 \
-soft_masking false \
-dust no \
-out "$wk_dir/blast/RE_blast.tsv"

#Extract the downstream of RE from all RE-contain contigs
python scripts/RE_blast_contig_fetch.py \
  -f "$wk_dir/blast/$sequencing_basename.fasta" \
  -t "$wk_dir/blast/RE_blast.tsv" \
  -s 200 \
  -o "$wk_dir/RE_contain_contigs.fasta"

echo "Mapping RE onto reference genome"
minimap2 -t 48 -ax map-ont "$genome_file" "$wk_dir/RE_contain_contigs.fasta" | samtools view -bS | samtools sort -o "$wk_dir/aligned_contigs.bam" - && samtools index "$wk_dir/aligned_contigs.bam"

# Map the raw reads onto the genome, then compute coverage
echo "Mapping raw fastq reads onto reference genome"
minimap2 -t 48 -ax map-ont "$genome_file" "$fastq_path" | samtools view -bS | samtools sort -o "$wk_dir/raw_reads_aligned.bam" - && samtools index "$wk_dir/raw_reads_aligned.bam"
samtools depth -o "$wk_dir/raw_reads_coverage.tsv" "$wk_dir/raw_reads_aligned.bam"
