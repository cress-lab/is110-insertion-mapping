#!/bin/bash

# Script to search an HMM against all of GTDB

# Set directories
GTDB_DIR="/shared/db/gtdb/latest/genomic_files_reps/protein_faa_reps"
ARCHAEA_DIR="$GTDB_DIR/archaea"
BACTERIA_DIR="$GTDB_DIR/bacteria"
HMM_FILE="$1"
OUTPUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over both domains
for DOMAIN in "$ARCHAEA_DIR" "$BACTERIA_DIR"; do
    for PROTEIN_FILE in "$DOMAIN"/*.faa; do
        # Get the genome name from the file (remove directory and extension)
        GENOME_NAME=$(basename "$PROTEIN_FILE" .faa)

        HMM_NAME=$(basename "$HMM_FILE" .hmm)

        # Define output file
        OUT_FILE="$OUTPUT_DIR/${GENOME_NAME}__${HMM_NAME}.domtblout"

        # Run hmmsearch with tblout
        hmmsearch --domtblout "$OUT_FILE" "$HMM_FILE" "$PROTEIN_FILE" > /dev/null
        done
    done
done
