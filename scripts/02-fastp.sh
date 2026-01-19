#!/bin/bash

# Load install fastp/0.23.4

SPECIES="Cannabis_sativa"
WORKDIR=$(dirname $0)
cd "$WORKDIR"

mkdir -p clean_reads fastpReports

# Loop through each line of metadata (skipping header if any)
while read -r RUN SAMPLE LAYOUT PLATFORM READTYPE EXTRA; do

    # Determine minimum length
    [[ "$READTYPE" == "short" ]] && LENGTH=30 || LENGTH=50

    # Locate raw reads (works for SRA / ENA / NGDC)
    RAW_READS=(RAW_DATA/"$RUN"*)
    
    if [[ "$LAYOUT" == "PE" ]]; then
        READ1="${RAW_READS[0]}"
        READ2="${RAW_READS[1]}"
        CLEAN1="clean_reads/${RUN}_1.fastq"
        CLEAN2="clean_reads/${RUN}_2.fastq"
    else
        READ1="${RAW_READS[0]}"
        CLEAN1="clean_reads/${RUN}.fastq"
    fi

    # Determine analysis type
    if   [[ "$LAYOUT" == "SE" && "$PLATFORM" == "Illumina" ]]; then TYPE=0
    elif [[ "$LAYOUT" == "PE" && "$PLATFORM" == "Illumina" ]]; then TYPE=1
    elif [[ "$LAYOUT" == "SE" && "$PLATFORM" == "DNB"      ]]; then TYPE=2
    elif [[ "$LAYOUT" == "PE" && "$PLATFORM" == "DNB"      ]]; then TYPE=3
    else
        echo "Error: Unknown layout/platform for $RUN"
        continue
    fi

    # Run fastp
    case $TYPE in
        0)  # SE Illumina
            fastp \
                -i "$READ1" \
                -o "$CLEAN1" \
                --qualified_quality_phred 20 \
                --unqualified_percent_limit 30 \
                --cut_front --cut_front_window_size 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 \
                --length_required $LENGTH \
                --json fastpReports/"$RUN".json
            ;;

        1)  # PE Illumina
            fastp \
                -i "$READ1" -I "$READ2" \
                -o "$CLEAN1" -O "$CLEAN2" \
                --qualified_quality_phred 20 \
                --unqualified_percent_limit 30 \
                --cut_front --cut_front_window_size 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 \
                --length_required $LENGTH \
                --json fastpReports/"$RUN".json
            ;;

        2)  # SE DNB
            fastp \
                -i "$READ1" \
                -o "$CLEAN1" \
                --adapter_fasta "$HOME/annotation_scripts/adaptors_DNBSEQ-T7.fasta" \
                --qualified_quality_phred 20 \
                --unqualified_percent_limit 30 \
                --cut_front --cut_front_window_size 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 \
                --length_required $LENGTH \
                --json fastpReports/"$RUN".json
            ;;

        3)  # PE DNB
            fastp \
                -i "$READ1" -I "$READ2" \
                -o "$CLEAN1" -O "$CLEAN2" \
                --adapter_fasta "$HOME/annotation_scripts/adaptors_DNBSEQ-T7.fasta" \
                --qualified_quality_phred 20 \
                --unqualified_percent_limit 30 \
                --cut_front --cut_front_window_size 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 \
                --length_required $LENGTH \
                --json fastpReports/"$RUN".json
            ;;
    esac

done < metadata/metadata.txt
