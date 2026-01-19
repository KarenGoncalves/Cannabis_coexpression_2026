#!/usr/bin/env bash

# Strandedness checker using RSeQC infer_experiment.py
# pip install how_are_we_stranded

set -euo pipefail

# ---------------------------------------------------------
# Project structure: scripts/ is inside the project folder
# ---------------------------------------------------------
PROJECT_DIR=$(cd $(dirname $0)/.. && pwd)

# Inputs
MAP_DIR="${PROJECT_DIR}/mapping_star"     # Directory with BAM files
GENOME_DIR="${PROJECT_DIR}/genome_index"  # Directory containing the GTF
GTF="${GENOME_DIR}/Cannabis_sativa.gtf"
BED="${GENOME_DIR}/Cannabis_sativa.bed"

# Output directory
OUT_DIR="${PROJECT_DIR}/strandness"
SUMMARY="${OUT_DIR}/strandness_summary.tsv"
mkdir -p "$OUT_DIR"

echo "RUNID    FRAC_FORWARD    FRAC_REVERSE    INTERPRETATION" > "$SUMMARY"

# ---------------------------------------------------------
# Create BED file if missing (required by RSeQC)
# ---------------------------------------------------------
if [[ ! -s "$BED" ]]; then
    echo "Creating BED from GTF..."
    gtf2bed.py < "$GTF" > "$BED"
fi

# ---------------------------------------------------------
# Process each BAM file sequentially
# ---------------------------------------------------------
for BAM in "${MAP_DIR}"/*.bam; do
    RUNID=$(basename "$BAM" .Aligned.sortedByCoord.out.bam)
    REPORT="${OUT_DIR}/${RUNID}.txt"

    echo "Processing $RUNID..."

    # Run RSeQC
    infer_experiment.py -r "$BED" -i "$BAM" > "$REPORT" 2>&1

    # Extract fractions
    FWD=$(grep '1++' "$REPORT" | awk -F': ' '{print $2}')
    REV=$(grep '1+-' "$REPORT" | awk -F': ' '{print $2}')

    # Simple interpretation
    if (( $(echo "$FWD > 0.8" | bc -l) )); then
        TYPE="Forward"
    elif (( $(echo "$REV > 0.8" | bc -l) )); then
        TYPE="Reverse"
    else
        TYPE="Unstranded/Unknown"
    fi

    # Add to summary
    echo -e "${RUNID}\t${FWD}\t${REV}\t${TYPE}" >> "$SUMMARY"
done

echo "Done. Summary: $SUMMARY"
