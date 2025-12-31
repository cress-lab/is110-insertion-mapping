#!/bin/bash

# Script to map sequences to 16S rRNA genes and then to genomes, 
# extracting unmapped reads at each step.
# USAGE: bash map_16S.sh <seqs_to_map.fasta> <16S_seqs.fasta> <genomes.fasta>

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

seqs_to_map="$1" # Fasta file of sequences to map
seqs_16s_fasta="$2" # Fasta file of 16S sequences to map against
genomes_fasta="$3" # Fasta file of genomes to map against

outdir=$(dirname "$seqs_to_map")/mapping_16s
mkdir -p "$outdir/blastdb"

mapped_to_16s_fasta="$outdir/16S_reads.fasta"
reads_16s_mapped_to_genomes_bam="$outdir/16S_mapped_to_genomes.bam"
mapped_to_genomes_bam="$outdir/non_16S_mapped_to_genomes.bam"
mapped_to_genomes_fasta="$outdir/$(basename "$mapped_to_genomes_bam" .bam).fasta"
unmapped_to_genomes_bam="$outdir/non_16S_mapped_to_genomes.unmapped.bam"
unmapped_to_genomes_fasta="$outdir/$(basename "$unmapped_to_genomes_bam" .bam).fasta"
blastn_out="$outdir/blastn_16s.txt"
# Use BLASTn to map the reads against the 16S database
makeblastdb -in $seqs_16s_fasta -dbtype nucl -out "$outdir/blastdb/16S"
blastn \
    -db "$outdir/blastdb/16S" \
    -query "$seqs_to_map" \
    -out "$blastn_out" \
    -outfmt 7
# python "$script_dir/process_16s_blastn.py" "$outdir" "$blastn_out" "$seqs_to_map"
python "$script_dir/process_16s_blastn.py" --blast_file "$blastn_out" --fasta_16s "$seqs_16s_fasta" --reads_file "$seqs_to_map"
unmapped_to_16s_fasta="$outdir/reads_no_16S_blastn_hit.fasta"
# Now, I’ll map these unmapped reads back to the genomes:
minimap2 -t 80 -ax map-ont $genomes_fasta $unmapped_to_16s_fasta |
samtools view -bS |
samtools sort -o $mapped_to_genomes_bam -
samtools index $mapped_to_genomes_bam
# Save the reads mapped to genomes to a fasta file
samtools fastq -F 2308 $mapped_to_genomes_bam | seqtk seq -A > "$mapped_to_genomes_fasta"
# Save the reads not mapped to genomes or 16S to a fasta file
samtools view -b -f 4 $mapped_to_genomes_bam > $unmapped_to_genomes_bam
samtools fastq $unmapped_to_genomes_bam | seqtk seq -A > "$unmapped_to_genomes_fasta"
# Finally, let's map the reads with 16S to the genome
minimap2 -t 80 -ax map-ont $genomes_fasta $mapped_to_16s_fasta |
samtools view -bS |
samtools sort -o $reads_16s_mapped_to_genomes_bam -
samtools index $reads_16s_mapped_to_genomes_bam
# Make a summary tsv of how many reads map where
N_16s="$(grep -c '>' $mapped_to_16s_fasta)"
N_off_target="$(grep -c '>' $mapped_to_genomes_fasta)"
N_unmapped="$(grep -c '>' $unmapped_to_genomes_fasta)"
echo -e "Category\tCount\n16S reads\t$N_16s\nNon-16S\t$N_off_target\nUnmapped to genome\t$N_unmapped" > "$outdir/16S_mapping_summary.tsv"
# Cleanup the blastdb
rm -r "$outdir/blastdb"
