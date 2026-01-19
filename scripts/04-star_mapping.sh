
#!/usr/bin/env bash
###########################################################################
# STAR READ MAPPING PIPELINE 
#
# - Reads metadata file and loops through each sample
# - Chooses appropriate STAR genome index based on read length
# - Handles SE/PE logic cleanly
# - Outputs BAM + STAR logs + stats
#
# REQUIREMENTS:
#   - genome_index_<RL>/ created by 03-star_index.sh
#   - cleaned reads inside:  clean_reads/
#   - STAR installed
###########################################################################

set -euo pipefail
THREADS=8

# -------------------------------------------------------------------------
# Project directory = parent folder of scripts/
# -------------------------------------------------------------------------
PROJECT_DIR=$(cd $(dirname $0)/.. && pwd)

# Input/output dirs
CLEAN_READS="${PROJECT_DIR}/clean_reads"
METADATA="${PROJECT_DIR}/metadata/metadata.txt"
OUTDIR="${PROJECT_DIR}/mapping_star"

mkdir -p ${OUTDIR}

# -------------------------------------------------------------------------
# Metadata columns expected:
# 1) RUN ID
# 2) Sample name
# 3) Layout                (SE|PE)
# 4) Platform              (Illumina|DNB)
# 5) Read length category  (short|long|specific length)
# 6) Actual read length    (optional numeric column)
#
# This script assumes column 6 contains a numeric read length.
# If not, modify accordingly.
# -------------------------------------------------------------------------

echo "=== STAR Mapping Pipeline ==="
echo "Project directory: $PROJECT_DIR"
echo "Metadata file:     $METADATA"
echo "Clean reads:       $CLEAN_READS"
echo "Output directory:  $OUTDIR"
echo

# -------------------------------------------------------------------------
# Main loop: read metadata line by line
# -------------------------------------------------------------------------
while read -r RUNID SAMPLENAME LAYOUT PLATFORM READTYPE READLEN; do
    # skip header or malformed lines
    [[ ${RUNID} == "Run" ]] && continue
    [[ -z ${RUNID} ]] && continue

    echo "---- Processing $RUNID ($LAYOUT, $PLATFORM, read length $READLEN) ----"

    # -------------------------
    # Identify clean reads
    # -------------------------
    if [[ ${LAYOUT} == "PE" ]]; then
        READ1="${CLEAN_READS}/${RUNID}_1.fastq"
        READ2="${CLEAN_READS}/${RUNID}_2.fastq"
    else
        READ1="${CLEAN_READS}/${RUNID}.fastq"
    fi

    # Check reads exist
    if [[ ! -s ${READ1} ]]; then
        echo "ERROR: Missing clean read file: $READ1"
        continue
    fi
    if [[ ${LAYOUT} == "PE" && ! -s ${READ2} ]]; then
        echo "ERROR: Missing clean read file: $READ2"
        continue
    fi

    # -------------------------
    # Determine which STAR index to use
    # Index directory naming: genome_index/Index_<readlength>
    # e.g. genome_index/Index_90
    # -------------------------
    STAR_INDEX="${PROJECT_DIR}/genome_index/Index_${READLEN}"

    if [[ ! -d ${STAR_INDEX} ]]; then
        echo "ERROR: No STAR index found for read length ${READLEN}"
        echo "Expected directory: $STAR_INDEX"
        echo "Make sure 03-star_index.sh created it."
        continue
    fi

    # -------------------------
    # Output files
    # -------------------------
    PREFIX="${OUTDIR}/${RUNID}"

    # -------------------------
    # Run STAR
    # -------------------------
    if [[ ${LAYOUT} == "PE" ]]; then
        STAR \
            --runThreadN ${THREADS} \
            --genomeDir ${STAR_INDEX} \
            --readFilesIn ${READ1} ${READ2} \
            --outFileNamePrefix "${PREFIX}." \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts
    else
        STAR \
            --runThreadN ${THREADS} \
            --genomeDir ${STAR_INDEX} \
            --readFilesIn ${READ1} \
            --outFileNamePrefix "${PREFIX}." \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts
    fi

    echo "Mapping completed for $RUNID"
    echo
done < ${METADATA}

echo "=== All samples processed ==="
