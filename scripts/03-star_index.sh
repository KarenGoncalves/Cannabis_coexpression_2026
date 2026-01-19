
#!/usr/bin/env bash
#####################################################
# Install star/2.7.11b
# Build STAR genome index for Cannabis sativa genome
#####################################################

##############################################
# Check STAR availability
##############################################
if ! command -v STAR >/dev/null 2>&1; then
  echo "Error: STAR not found in PATH. Please install STAR." >&2
  exit 1
fi

set -euo pipefail
THREADS=8

# Project root = parent folder of scripts/
PROJECT_DIR=$(cd $(dirname $0)/.. && pwd)

# Paths
GENOME_DIR="${PROJECT_DIR}/genome_index"     # base directory (also where FASTA/GTF are stored)
SPECIES="Cannabis_sativa"
FASTA="${GENOME_DIR}/${SPECIES}.fa"          # target uncompressed FASTA
GTF="${GENOME_DIR}/${SPECIES}.gtf"           # target uncompressed GTF

# Read lengths you want to support; sjdbOverhang = read_length - 1
READ_LENGTHS=(30 70 90 144)

# Source URLs (Ensembl Plants)
FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/cannabis_sativa_female/dna/Cannabis_sativa_female.cs10.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/cannabis_sativa_female/Cannabis_sativa_female.cs10.59.gtf.gz"

mkdir -p ${GENOME_DIR}

##############################################
# Helper: Download + Decompress (.gz -> plain)
# Usage: download_and_decompress URL OUT_PLAIN
# Behavior:
#   - If OUT_PLAIN exists and non-empty > skip
#   - Else download .gz to OUT_PLAIN.gz (if missing)
#   - Decompress into OUT_PLAIN (atomic via temp file)
##############################################
download_and_decompress() {
  local url=${1}
  local out_plain=${2}
  local out_gz="${out_plain}.gz"

  # If the final plain file is already present, skip
  if [[ -s ${out_plain} ]]; then
    echo "Present (plain): $out_plain"
    return
  fi

  # Download the gz if missing
  if [[ ! -s ${out_gz} ]]; then
    echo "Downloading: $url"
    if command -v wget >/dev/null 2>&1; then
      wget -q ${url} -O ${out_gz}
    elif command -v curl >/dev/null 2>&1; then
      curl -sSL ${url} -o ${out_gz}
    else
      echo "Error: Neither wget nor curl found." >&2
      exit 1
    fi
  else
    echo "Present (gz): $out_gz"
  fi

  # Decompress to a temp file then move atomically
  echo "Decompressing: $out_gz → $out_plain"
  local tmp
  tmp=$(mktemp ${out_plain}.XXXXXX)
  gzip -dc ${out_gz} > ${tmp}
  mv -f ${tmp} ${out_plain}
}

##############################################
# Fetch inputs (download + decompress)
##############################################
download_and_decompress ${FASTA_URL} ${FASTA}
download_and_decompress ${GTF_URL}   ${GTF}


##############################################
# Build STAR indices (one per read length)
##############################################
echo "Building STAR indices for ${SPECIES}..."
echo "Threads     : ${THREADS}"
echo "FASTA       : ${FASTA}"
echo "GTF         : ${GTF}"
echo "Base out    : ${GENOME_DIR}"
echo "Read lengths: ${READ_LENGTHS[@]}"
echo

for RL in "${READ_LENGTHS[@]}"; do
  # Validate read length
  if ! [[ ${RL} =~ ^[0-9]+$ ]] || (( RL <= 0 )); then
    echo "Skipping invalid read length: $RL" >&2
    continue
  fi

  OVERHANG=$((RL - 1))
  OUTDIR="${GENOME_DIR}_${RL}"

  echo "==> Building STAR index: read length ${RL} (sjdbOverhang=${OVERHANG})"
  echo "    Output dir: ${OUTDIR}"
  mkdir -p ${OUTDIR}

  STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${OUTDIR} \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang ${OVERHANG}

  echo "Done: ${OUTDIR}"
  echo
done

echo "All requested STAR indices built."
