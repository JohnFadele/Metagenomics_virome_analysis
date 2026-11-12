#!/usr/bin/env bash

###############################################################################
#                 GENERIC REFERENCE-BASED ASSEMBLY PIPELINE
###############################################################################
#
# PURPOSE:
#
#   1. Detect all paired FASTQ files in INPUT_DIR
#   2. Support .fastq and .fastq.gz
#   3. Support arbitrary suffixes after _R1 / _R2
#   4. Require exact R1/R2 filename matching except R1 vs R2
#   5. Detect all reference FASTA files in REF_DIR
#   6. Index references where necessary
#   7. Map every sample against every reference
#   8. Generate consensus sequences using iVar
#   9. Calculate resolved bases
#  10. Calculate percentage of reference resolved
#  11. Organise all consensus FASTA files
#  12. Produce one wide-format CSV
#  13. Remove unnecessary intermediate BAM files
#  14. Resume previously completed sample/reference analyses
#
# IMPORTANT:
#
#   - INPUT FASTQ FILES ARE NEVER MODIFIED.
#   - INPUT FASTQ FILES ARE NEVER RENAMED.
#   - INPUT FASTQ FILES ARE NEVER DELETED.
#   - Existing completed assemblies are skipped.
#   - Incomplete assemblies are rerun.
#
###############################################################################

set -uo pipefail

IFS=$'\n\t'


###############################################################################
#                         USER SETTINGS
###############################################################################

# Number of CPU threads

THREADS=8


# ---------------------------------------------------------------------------
# INPUT FASTQ DIRECTORY
# ---------------------------------------------------------------------------

INPUT_DIR="/home/fadele/ppr_extra/results/ppr_reads"


# ---------------------------------------------------------------------------
# REFERENCE DIRECTORY
#
# The script automatically detects:
#
# *.fa
# *.fasta
# *.fna
#
# ---------------------------------------------------------------------------

REF_DIR="/home/fadele/references/PPRV"


# ---------------------------------------------------------------------------
# MAIN OUTPUT DIRECTORY
# ---------------------------------------------------------------------------

OUT_DIR="/home/fadele/ppr_extra/results/ppr_assembly"


###############################################################################
#                         OUTPUT DIRECTORIES
###############################################################################

ASSEMBLY_DIR="$OUT_DIR/assembly"

ALL_FASTA_DIR="$OUT_DIR/all_fasta"

LOG_DIR="$OUT_DIR/logs"

SUMMARY="$OUT_DIR/assembly_summary_wide.csv"


###############################################################################
#                         CREATE DIRECTORIES
###############################################################################

mkdir -p "$OUT_DIR"
mkdir -p "$ASSEMBLY_DIR"
mkdir -p "$ALL_FASTA_DIR"
mkdir -p "$LOG_DIR"


###############################################################################
#                         CHECK REQUIRED PROGRAMS
###############################################################################

echo
echo "============================================================"
echo "       CHECKING REQUIRED PROGRAMS"
echo "============================================================"

for PROGRAM in bwa samtools ivar
do

    if ! command -v "$PROGRAM" >/dev/null 2>&1
    then

        echo
        echo "[ERROR] $PROGRAM was not found."
        echo
        echo "Please activate the appropriate conda environment."
        echo

        exit 1

    fi

    echo "[OK] $PROGRAM"

done


###############################################################################
#                         CHECK INPUT DIRECTORY
###############################################################################

if [[ ! -d "$INPUT_DIR" ]]
then

    echo
    echo "[ERROR] Input FASTQ directory does not exist:"
    echo "$INPUT_DIR"
    echo

    exit 1

fi


###############################################################################
#                         CHECK REFERENCE DIRECTORY
###############################################################################

if [[ ! -d "$REF_DIR" ]]
then

    echo
    echo "[ERROR] Reference directory does not exist:"
    echo "$REF_DIR"
    echo

    exit 1

fi


###############################################################################
#                         DETECT REFERENCES
###############################################################################

echo
echo "============================================================"
echo "       DETECTING REFERENCE FASTA FILES"
echo "============================================================"

REFERENCES=()

while IFS= read -r REF
do

    REFERENCES+=("$REF")

done < <(

    find "$REF_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name "*.fa" \
            -o -name "*.fasta" \
            -o -name "*.fna" \
        \) |
    sort

)


###############################################################################
#                         CHECK REFERENCES
###############################################################################

if [[ ${#REFERENCES[@]} -eq 0 ]]
then

    echo
    echo "[ERROR] No reference FASTA files were found."
    echo
    echo "Reference directory:"
    echo "$REF_DIR"
    echo

    exit 1

fi


echo
echo "Detected ${#REFERENCES[@]} reference FASTA files:"
echo

for REF in "${REFERENCES[@]}"
do

    echo "  $(basename "$REF")"

done


###############################################################################
#                         INDEX REFERENCES
###############################################################################

echo
echo "============================================================"
echo "       CHECKING / BUILDING REFERENCE INDEXES"
echo "============================================================"


for REF in "${REFERENCES[@]}"
do

    REF_NAME=$(basename "$REF")

    echo
    echo "Reference:"
    echo "$REF_NAME"


    ###########################################################################
    # BWA INDEX
    ###########################################################################

    if [[ ! -f "${REF}.bwt" ]]
    then

        echo "[INFO] BWA index not found."
        echo "[INFO] Building BWA index..."

        if ! bwa index "$REF"
        then

            echo
            echo "[ERROR] Failed to index:"
            echo "$REF"
            echo

            exit 1

        fi

        echo "[OK] BWA index created."

    else

        echo "[OK] BWA index already exists."

    fi


    ###########################################################################
    # SAMTOOLS FASTA INDEX
    ###########################################################################

    if [[ ! -f "${REF}.fai" ]]
    then

        echo "[INFO] FASTA index not found."
        echo "[INFO] Building FASTA index..."

        if ! samtools faidx "$REF"
        then

            echo
            echo "[ERROR] Failed to create FASTA index:"
            echo "$REF"
            echo

            exit 1

        fi

        echo "[OK] FASTA index created."

    else

        echo "[OK] FASTA index already exists."

    fi

done


###############################################################################
#                         CREATE CSV HEADER
###############################################################################

echo
echo "============================================================"
echo "       CREATING SUMMARY FILE"
echo "============================================================"


HEADER="S/N,Sample_ID"


for REF in "${REFERENCES[@]}"
do

    REF_NAME=$(basename "$REF")

    REF_ID="$REF_NAME"

    REF_ID="${REF_ID%.fasta}"
    REF_ID="${REF_ID%.fa}"
    REF_ID="${REF_ID%.fna}"


    HEADER="${HEADER},${REF_ID}_Resolved_Bases"

    HEADER="${HEADER},${REF_ID}_Ref_Len_bp"

    HEADER="${HEADER},${REF_ID}_Pct_Resolved"

done


echo "$HEADER" > "$SUMMARY"


###############################################################################
#                         DETECT SAMPLES
###############################################################################
#
# Supported examples:
#
#   SAMPLE_R1.fastq
#   SAMPLE_R2.fastq
#
#   SAMPLE_R1.fastq.gz
#   SAMPLE_R2.fastq.gz
#
#   SAMPLE_R1.trimmed.fastq
#   SAMPLE_R2.trimmed.fastq
#
#   SAMPLE_R1.trimmed.fastq.gz
#   SAMPLE_R2.trimmed.fastq.gz
#
#   SAMPLE_R1.clean.fastq
#   SAMPLE_R2.clean.fastq
#
#   SAMPLE_R1.dehosted.fastq.gz
#   SAMPLE_R2.dehosted.fastq.gz
#
#   and arbitrary suffixes following _R1 / _R2.
#
###############################################################################

echo
echo "============================================================"
echo "       DETECTING PAIRED FASTQ FILES"
echo "============================================================"

R1_FILES=()

while IFS= read -r R1
do

    R1_FILES+=("$R1")

done < <(

    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -name "*_R1.fastq" \
            -o -name "*_R1.fastq.gz" \
            -o -name "*_R1.*.fastq" \
            -o -name "*_R1.*.fastq.gz" \
        \) |
    sort

)


###############################################################################
#                         CHECK SAMPLES
###############################################################################

if [[ ${#R1_FILES[@]} -eq 0 ]]
then

    echo
    echo "[ERROR] No paired-end R1 FASTQ files were found."
    echo
    echo "Input directory:"
    echo "$INPUT_DIR"
    echo

    exit 1

fi


echo
echo "Detected ${#R1_FILES[@]} R1 FASTQ files."
echo


###############################################################################
#                         FUNCTION:
#                         REFERENCE LENGTH
###############################################################################

get_reference_length()
{

    local REF="$1"

    awk '

        !/^>/ {
            total += length($0)
        }

        END {
            print total + 0
        }

    ' "$REF"

}


###############################################################################
#                         FUNCTION:
#                         RESOLVED BASES
###############################################################################

get_resolved_bases()
{

    local CONSENSUS="$1"


    if [[ ! -s "$CONSENSUS" ]]
    then

        echo "0"

        return

    fi


    grep -v "^>" "$CONSENSUS" |
        tr -cd 'ACGTacgt' |
        wc -c

}


###############################################################################
#                         FUNCTION:
#                         CALCULATE PERCENT
###############################################################################

calculate_percent()
{

    local RESOLVED="$1"
    local REF_LEN="$2"


    if [[ "$REF_LEN" -gt 0 ]]
    then

        awk \
            -v resolved="$RESOLVED" \
            -v reflength="$REF_LEN" \
            'BEGIN {
                printf "%.2f", (resolved / reflength) * 100
            }'

    else

        echo "0.00"

    fi

}


###############################################################################
#                         SAMPLE LOOP
###############################################################################

SN=1


for R1 in "${R1_FILES[@]}"
do

    ###########################################################################
    # R1 FILENAME
    ###########################################################################

    R1_BASE=$(basename "$R1")


    ###########################################################################
    # CREATE EXACT R2 FILENAME
    #
    # Everything remains identical except _R1 becomes _R2.
    ###########################################################################

    R2_BASE="${R1_BASE/_R1/_R2}"

    R2="$INPUT_DIR/$R2_BASE"


    ###########################################################################
    # SAMPLE NAME
    #
    # Everything before the first _R1 is treated as the sample identifier.
    ###########################################################################

    SAMPLE="${R1_BASE%%_R1*}"


    echo
    echo
    echo "============================================================"
    echo "PROCESSING SAMPLE: $SAMPLE"
    echo "============================================================"

    echo
    echo "[INPUT]"
    echo "R1: $R1"
    echo "R2: $R2"


    ###########################################################################
    # CHECK R2
    ###########################################################################

    if [[ ! -f "$R2" ]]
    then

        echo
        echo "[WARNING] Matching R2 file was not found."

        echo "[WARNING] Expected:"
        echo "$R2"

        echo "[WARNING] Skipping sample."

        continue

    fi


    ###########################################################################
    # SAMPLE OUTPUT DIRECTORY
    ###########################################################################

    SAMPLE_OUT="$ASSEMBLY_DIR/$SAMPLE"

    mkdir -p "$SAMPLE_OUT"


    ###########################################################################
    # SAMPLE LOG
    ###########################################################################

    SAMPLE_LOG="$LOG_DIR/${SAMPLE}.log"


    echo "Sample: $SAMPLE" > "$SAMPLE_LOG"
    echo "R1: $R1" >> "$SAMPLE_LOG"
    echo "R2: $R2" >> "$SAMPLE_LOG"
    echo "Started: $(date)" >> "$SAMPLE_LOG"


    ###########################################################################
    # START CSV ROW
    ###########################################################################

    ROW="$SN,$SAMPLE"


    ###########################################################################
    # LOOP THROUGH EVERY REFERENCE
    ###########################################################################

    for REF in "${REFERENCES[@]}"
    do

        REF_NAME=$(basename "$REF")


        #######################################################################
        # REFERENCE ID
        #######################################################################

        REF_ID="$REF_NAME"

        REF_ID="${REF_ID%.fasta}"
        REF_ID="${REF_ID%.fa}"
        REF_ID="${REF_ID%.fna}"


        echo
        echo "------------------------------------------------------------"
        echo "Sample    : $SAMPLE"
        echo "Reference : $REF_ID"
        echo "------------------------------------------------------------"


        #######################################################################
        # REFERENCE LENGTH
        #######################################################################

        REF_LEN=$(get_reference_length "$REF")

        echo "[INFO] Reference length: $REF_LEN bp"


        #######################################################################
        # OUTPUT FILES
        #######################################################################

        BAM="$SAMPLE_OUT/${SAMPLE}_${REF_ID}_sorted.bam"

        BAI="${BAM}.bai"

        CONS_PREFIX="$SAMPLE_OUT/${SAMPLE}_${REF_ID}"

        CONSENSUS="$CONS_PREFIX.fa"

        CENTRAL_FASTA="$ALL_FASTA_DIR/${SAMPLE}_${REF_ID}.fa"


        #######################################################################
        # RESUME CHECK
        #
        # A non-empty consensus FASTA is treated as a successfully completed
        # sample/reference analysis.
        #######################################################################

        if [[ -s "$CONSENSUS" ]]
        then

            echo "[SKIP] Completed consensus already exists."

            echo "[SKIP] Consensus:"
            echo "$CONSENSUS"


            ###################################################################
            # RECOVER SUMMARY VALUES FROM EXISTING CONSENSUS
            ###################################################################

            RESOLVED=$(get_resolved_bases "$CONSENSUS")

            PERCENT=$(calculate_percent "$RESOLVED" "$REF_LEN")


            echo
            echo "[EXISTING RESULT]"
            echo "Reference length : $REF_LEN bp"
            echo "Resolved bases   : $RESOLVED bp"
            echo "Percent resolved : $PERCENT %"


            ###################################################################
            # ENSURE CENTRAL FASTA COPY EXISTS
            ###################################################################

            if [[ ! -s "$CENTRAL_FASTA" ]]
            then

                echo "[INFO] Central FASTA copy missing."
                echo "[INFO] Creating central FASTA copy."

                cp "$CONSENSUS" "$CENTRAL_FASTA"

            else

                echo "[OK] Central FASTA copy already exists."

            fi


            ###################################################################
            # ADD RESULTS TO CSV ROW
            ###################################################################

            ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


            ###################################################################
            # MOVE TO NEXT REFERENCE
            ###################################################################

            continue

        fi


        #######################################################################
        # REMOVE STALE INTERMEDIATE FILES
        #
        # If a consensus does not exist, any incomplete BAM/BAI from a previous
        # failed run is removed and the reference is rebuilt.
        #######################################################################

        if [[ -f "$BAM" || -f "$BAI" ]]
        then

            echo "[INFO] Incomplete intermediate files detected."

            echo "[INFO] Removing stale BAM/BAI before rerunning."

            rm -f "$BAM"
            rm -f "$BAI"

        fi


        #######################################################################
        # MAP READS
        #######################################################################

        echo "[INFO] Mapping reads..."

        echo "[INFO] R1: $R1"
        echo "[INFO] R2: $R2"


        #######################################################################
        # BWA MEM
        #
        # BWA can read both uncompressed FASTQ and gzip-compressed FASTQ
        # directly, so no decompression is performed.
        #######################################################################

        if ! bwa mem \
            -t "$THREADS" \
            "$REF" \
            "$R1" \
            "$R2" \
            2>> "$SAMPLE_LOG" |

            samtools view \
                -b \
                -F 4 \
                - |

            samtools sort \
                -o "$BAM" \
                -
        then

            echo "[WARNING] Mapping failed."
            echo "[WARNING] Reference: $REF_ID"

            RESOLVED=0

            PERCENT="0.00"


            ###################################################################
            # ADD FAILURE VALUES TO SAMPLE ROW
            ###################################################################

            ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


            ###################################################################
            # REMOVE FAILED BAM
            ###################################################################

            rm -f "$BAM"
            rm -f "$BAI"


            continue

        fi


        #######################################################################
        # INDEX BAM
        #######################################################################

        echo "[INFO] Indexing BAM..."


        if ! samtools index "$BAM"
        then

            echo "[WARNING] BAM indexing failed."


            RESOLVED=0

            PERCENT="0.00"


            ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


            rm -f "$BAM"
            rm -f "$BAI"


            continue

        fi


        #######################################################################
        # GENERATE CONSENSUS
        #######################################################################

        echo "[INFO] Generating iVar consensus..."


        if ! samtools mpileup \
            -aa \
            -A \
            -d 0 \
            -Q 0 \
            "$BAM" |

            ivar consensus \
                -p "$CONS_PREFIX" \
                -t 0.5 \
                -m 10
        then

            echo "[WARNING] Consensus generation failed."


            RESOLVED=0

            PERCENT="0.00"


            ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


            rm -f "$BAM"
            rm -f "$BAI"


            continue

        fi


        #######################################################################
        # CHECK CONSENSUS
        #######################################################################

        if [[ ! -s "$CONSENSUS" ]]
        then

            echo "[WARNING] Consensus FASTA was not produced."
            echo "[WARNING] Reference: $REF_ID"


            RESOLVED=0

            PERCENT="0.00"


            ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


            rm -f "$BAM"
            rm -f "$BAI"


            continue

        fi


        #######################################################################
        # COUNT RESOLVED BASES
        #######################################################################

        RESOLVED=$(get_resolved_bases "$CONSENSUS")


        #######################################################################
        # CALCULATE PERCENTAGE RESOLVED
        #######################################################################

        PERCENT=$(calculate_percent "$RESOLVED" "$REF_LEN")


        #######################################################################
        # DISPLAY RESULTS
        #######################################################################

        echo
        echo "[RESULT]"
        echo "Reference length : $REF_LEN bp"
        echo "Resolved bases   : $RESOLVED bp"
        echo "Percent resolved : $PERCENT %"


        #######################################################################
        # ADD RESULTS TO SAME SAMPLE ROW
        #######################################################################

        ROW="${ROW},${RESOLVED},${REF_LEN},${PERCENT}"


        #######################################################################
        # COPY CONSENSUS TO CENTRAL FASTA DIRECTORY
        #######################################################################

        cp "$CONSENSUS" \
            "$CENTRAL_FASTA"


        #######################################################################
        # REMOVE INTERMEDIATE BAM FILES
        #######################################################################

        rm -f "$BAM"

        rm -f "$BAI"


        echo "[OK] Reference completed: $REF_ID"

    done


    ###########################################################################
    # WRITE COMPLETE SAMPLE ROW
    ###########################################################################

    echo "$ROW" >> "$SUMMARY"


    ###########################################################################
    # SAMPLE LOG
    ###########################################################################

    echo "Finished: $(date)" >> "$SAMPLE_LOG"


    ###########################################################################
    # SAMPLE NUMBER
    ###########################################################################

    SN=$((SN + 1))


    echo
    echo "[DONE] Sample completed: $SAMPLE"

done


###############################################################################
#                         FINAL REPORT
###############################################################################

echo
echo
echo "============================================================"
echo "             ASSEMBLY PIPELINE COMPLETE"
echo "============================================================"

echo
echo "Input FASTQ directory:"
echo "$INPUT_DIR"

echo
echo "Reference directory:"
echo "$REF_DIR"

echo
echo "Assembly directory:"
echo "$ASSEMBLY_DIR"

echo
echo "All FASTA files:"
echo "$ALL_FASTA_DIR"

echo
echo "Assembly summary:"
echo "$SUMMARY"


###############################################################################
#                         SUMMARY COUNTS
###############################################################################

DATA_ROWS=$(tail -n +2 "$SUMMARY" | wc -l)


echo
echo "Number of samples processed: $DATA_ROWS"

echo
echo "============================================================"
echo "                  PIPELINE FINISHED"
echo "============================================================"

echo
echo "Finished: $(date)"
