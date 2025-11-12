#!/bin/bash

# Script to map sequences to 16S rRNA genes and then to genomes, 
# extracting unmapped reads at each step.
# USAGE: bash map_16S.sh <seqs_to_map.fasta> <16S_seqs.fasta> <genomes.fasta>

seqs_to_map="$1" # Fasta file of sequences to map
seqs_16s_fasta="$2" # Fasta file of 16S sequences to map against
genomes_fasta="$3" # Fasta file of genomes to map against

outdir=$(dirname "$seqs_to_map")/mapping_16s
mkdir -p "$outdir"

mapped_to_16s_bam="$outdir/16S_map.bam"
unmapped_to_16s_bam="$outdir/16S_map.unmapped.bam"
unmapped_to_16s_fasta="$outdir/$(basename "$unmapped_to_16s_bam" .bam).fasta"
mapped_to_genomes_bam="$outdir/non_16S_mapped_to_genomes.bam"
unmapped_to_genomes_bam="$outdir/non_16S_mapped_to_genomes.unmapped.bam"
# mapped_to_genomes_fasta=$(basename "$mapped_to_genomes_bam" .bam).fasta

# Do the mapping
minimap2 -t 80 -ax map-ont $seqs_16s_fasta $seqs_to_map |
samtools view -bS |
samtools sort -o $mapped_to_16s_bam -
samtools index $mapped_to_16s_bam
# Extract the reads mapped to 16S to a fasta file
# 4 → unmapped
# 256 → secondary
# 2048 → supplementary
# 4 + 256 + 2048 = 2308
samtools fastq -F 2308 $mapped_to_16s_bam | seqtk seq -A > "$outdir/$(basename "$mapped_to_16s_bam" .bam).fasta"
# Extract the unmapped reads to a fasta file
samtools view -b -f 4 $mapped_to_16s_bam > $unmapped_to_16s_bam
samtools fastq $unmapped_to_16s_bam | seqtk seq -A > $unmapped_to_16s_fasta
# Now, I’ll map these reads back to the genomes:
minimap2 -t 80 -ax map-ont $genomes_fasta $unmapped_to_16s_fasta |
samtools view -bS |
samtools sort -o $mapped_to_genomes_bam -
samtools index $mapped_to_genomes_bam
# Save the reads mapped to genomes to a fasta file
samtools fastq -F 2308 $mapped_to_genomes_bam | seqtk seq -A > "$outdir/$(basename "$mapped_to_genomes_bam" .bam).fasta"
# Save the reads not mapped to genomes or 16S to a fasta file
samtools view -b -f 4 $mapped_to_genomes_bam > $unmapped_to_genomes_bam
samtools fastq $unmapped_to_genomes_bam | seqtk seq -A > "$outdir/$(basename "$unmapped_to_genomes_bam" .bam).fasta"