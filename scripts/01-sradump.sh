#!/bin/sh

# Install  sra-toolkit/3.0.0
NCPUS=8

# Extract accession list from metadata (first column)
acc=($(cut -f 1  metadata/metadata.txt))

# Ensure output directory exists
mkdir -p RAW_DATA/

# Loop through each accession and download the reads
for accession in "${acc[@]}"; do
    fasterq-dump \
        --seq-defline '@$ac.$sn./$ri' \   # Define FASTQ header format
        --split-files \                   # Produce _1.fastq and _2.fastq for paired-end reads
        --outdir RAW_DATA/ \              # Output directory
        -e $NCPUS \                       # Use $NCPUS CPU threads
        "$accession" &                    # Run in background; process next sample
done

# OPTIONAL: wait for all background jobs before exiting
wait
echo "All downloads completed."
