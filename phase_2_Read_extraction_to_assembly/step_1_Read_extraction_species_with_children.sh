#!/usr/bin/env bash

###############################################################################
# TAXON + ALL DESCENDANTS READ EXTRACTION
#
# Extracts reads assigned to a specified taxon ID and ALL of its descendants
# from paired-end FASTQ files.
#
# PAIRING RULE:
#   R1 and R2 must be identical apart from the R1/R2 designation.
#
# RESUME / SKIP:
#   - Both output R1 and R2 exist and are non-empty -> SKIP
#   - Only one output exists -> remove incomplete output and rerun
#   - Neither output exists -> run extraction
#
# INPUT FILES ARE NEVER RENAMED, MOVED, DELETED, OR MODIFIED.
###############################################################################

set -u
IFS=$'\n\t'

###############################################################################
# USER SETTINGS
###############################################################################

# KrakenTools
KRAKENTOOLS="/home/fadele/KrakenTools/extract_kraken_reads.py"

# Input FASTQ directory
INPUT_DIR="/home/fadele/ppr_extra/results/dehosted"

# Kraken2 output directory
KRAKEN_DIR="/home/fadele/ppr_extra/results/kraken2"

# Output directory for extracted reads
OUTPUT_DIR="/home/fadele/ppr_extra/results/ppr_reads"

###############################################################################
# TARGET TAXON
###############################################################################

# Official taxonomic identifier used for extraction
TAXON_ID=3052343

# Human-readable name of the target taxon
TAXON_NAME="ppr"

###############################################################################
# CREATE OUTPUT DIRECTORY
###############################################################################

mkdir -p "$OUTPUT_DIR"

###############################################################################
# SUMMARY FILE
###############################################################################

SUMMARY="$OUTPUT_DIR/extraction_summary.csv"

echo "Sample,Status,R1_reads,R2_reads" > "$SUMMARY"

###############################################################################
# START
###############################################################################

echo
echo "============================================================"
echo "       TAXON + DESCENDANTS READ EXTRACTION"
echo "============================================================"
echo
echo "Target taxid : $TAXON_ID"
echo "Target name  : $TAXON_NAME"
echo
echo "The extraction includes:"
echo "  - Target taxon"
echo "  - All children"
echo "  - All grandchildren"
echo "  - All descendants at every level"
echo
echo "Input directory:"
echo "$INPUT_DIR"
echo
echo "Kraken directory:"
echo "$KRAKEN_DIR"
echo
echo "Output directory:"
echo "$OUTPUT_DIR"
echo
echo "============================================================"

###############################################################################
# COUNTERS
###############################################################################

FOUND=0
PROCESSED=0
SKIPPED=0
NO_READS=0
FAILED=0
MISSING_R2=0

###############################################################################
# FIND ALL R1 FILES
#
# Accepts:
#
#   *_R1.fastq
#   *_R1.fastq.gz
#   *_R1.<anything>.fastq
#   *_R1.<anything>.fastq.gz
#
###############################################################################

for R1 in "$INPUT_DIR"/*_R1.fastq \
           "$INPUT_DIR"/*_R1.fastq.gz \
           "$INPUT_DIR"/*_R1.*.fastq \
           "$INPUT_DIR"/*_R1.*.fastq.gz
do

    ###########################################################################
    # IGNORE UNMATCHED GLOB PATTERNS
    ###########################################################################

    [[ -f "$R1" ]] || continue

    FOUND=$((FOUND + 1))

    ###########################################################################
    # GET R1 FILENAME
    ###########################################################################

    R1_BASENAME=$(basename "$R1")

    ###########################################################################
    # IDENTIFY EXACT MATCHING R2
    #
    # Everything except R1/R2 must remain identical.
    ###########################################################################

    R2_BASENAME="${R1_BASENAME/_R1/_R2}"

    R2="$INPUT_DIR/$R2_BASENAME"

    ###########################################################################
    # SAMPLE NAME
    #
    # Everything before _R1 is treated as the sample name.
    ###########################################################################

    SAMPLE="${R1_BASENAME%%_R1*}"

    ###########################################################################
    # DISPLAY SAMPLE
    ###########################################################################

    echo
    echo "============================================================"
    echo "Checking sample: $SAMPLE"
    echo "============================================================"

    echo "R1:"
    echo "$R1"

    echo
    echo "Expected R2:"
    echo "$R2"

    ###########################################################################
    # CHECK R2
    ###########################################################################

    if [[ ! -f "$R2" ]]
    then

        echo
        echo "[WARNING] Matching R2 file not found."
        echo "[WARNING] Skipping this sample."

        echo "$SAMPLE,MISSING_R2,0,0" >> "$SUMMARY"

        MISSING_R2=$((MISSING_R2 + 1))

        continue

    fi

    ###########################################################################
    # CHECK COMPRESSION STATUS
    #
    # R1 and R2 must both be:
    #
    #   .fastq
    #
    # OR
    #
    #   .fastq.gz
    ###########################################################################

    R1_GZ=0
    R2_GZ=0

    [[ "$R1" == *.gz ]] && R1_GZ=1
    [[ "$R2" == *.gz ]] && R2_GZ=1

    if [[ "$R1_GZ" -ne "$R2_GZ" ]]
    then

        echo
        echo "[WARNING] R1 and R2 have different compression status."
        echo "[WARNING] R1: $R1"
        echo "[WARNING] R2: $R2"
        echo "[WARNING] Skipping sample."

        echo "$SAMPLE,COMPRESSION_MISMATCH,0,0" >> "$SUMMARY"

        FAILED=$((FAILED + 1))

        continue

    fi

    ###########################################################################
    # KRAKEN2 FILES
    #
    # Actual Kraken2 naming convention:
    #
    #   SAMPLE.kraken.output
    #   SAMPLE.report.txt
    ###########################################################################

    KRAKEN_OUTPUT="$KRAKEN_DIR/${SAMPLE}.kraken.output"

    KRAKEN_REPORT="$KRAKEN_DIR/${SAMPLE}.report.txt"

    ###########################################################################
    # OUTPUT FASTQ FILES
    ###########################################################################

    OUTPUT_R1="$OUTPUT_DIR/${SAMPLE}_R1.fastq"

    OUTPUT_R2="$OUTPUT_DIR/${SAMPLE}_R2.fastq"

    ###########################################################################
    # RESUME / SKIP CHECK
    #
    # BOTH files must exist AND contain data before the sample is skipped.
    ###########################################################################

    if [[ -s "$OUTPUT_R1" && -s "$OUTPUT_R2" ]]
    then

        echo
        echo "[SKIP] $SAMPLE has already been successfully extracted."
        echo "       Existing R1: $OUTPUT_R1"
        echo "       Existing R2: $OUTPUT_R2"

        R1_READS=$(awk 'END {print NR/4}' "$OUTPUT_R1")
        R2_READS=$(awk 'END {print NR/4}' "$OUTPUT_R2")

        echo "       R1 reads: $R1_READS"
        echo "       R2 reads: $R2_READS"

        echo "$SAMPLE,SKIPPED,$R1_READS,$R2_READS" >> "$SUMMARY"

        SKIPPED=$((SKIPPED + 1))

        continue

    fi

    ###########################################################################
    # HANDLE INCOMPLETE PREVIOUS EXTRACTION
    #
    # If only one output exists, or either output is empty, remove the
    # incomplete output pair before rerunning.
    ###########################################################################

    if [[ -e "$OUTPUT_R1" || -e "$OUTPUT_R2" ]]
    then

        echo
        echo "[WARNING] Incomplete previous extraction detected."
        echo "[INFO] Removing incomplete output files before rerunning."

        rm -f "$OUTPUT_R1" "$OUTPUT_R2"

    fi

    ###########################################################################
    # CHECK KRAKEN OUTPUT
    ###########################################################################

    if [[ ! -f "$KRAKEN_OUTPUT" ]]
    then

        echo
        echo "[WARNING] Kraken output not found:"
        echo "$KRAKEN_OUTPUT"

        echo "$SAMPLE,NO_KRAKEN_OUTPUT,0,0" >> "$SUMMARY"

        FAILED=$((FAILED + 1))

        continue

    fi

    ###########################################################################
    # CHECK KRAKEN REPORT
    ###########################################################################

    if [[ ! -f "$KRAKEN_REPORT" ]]
    then

        echo
        echo "[WARNING] Kraken report not found:"
        echo "$KRAKEN_REPORT"

        echo "$SAMPLE,NO_KRAKEN_REPORT,0,0" >> "$SUMMARY"

        FAILED=$((FAILED + 1))

        continue

    fi

    ###########################################################################
    # EXTRACT TARGET TAXON + ALL DESCENDANTS
    ###########################################################################

    echo
    echo "[INFO] Extracting:"
    echo "       Taxid $TAXON_ID"
    echo "       $TAXON_NAME"
    echo "       + ALL descendants"
    echo

    if ! python3 "$KRAKENTOOLS" \
        -k "$KRAKEN_OUTPUT" \
        --report "$KRAKEN_REPORT" \
        -s "$R1" \
        -s2 "$R2" \
        -t "$TAXON_ID" \
        --include-children \
        --fastq-output \
        -o "$OUTPUT_R1" \
        -o2 "$OUTPUT_R2"
    then

        echo
        echo "[WARNING] Extraction failed for $SAMPLE."

        rm -f "$OUTPUT_R1" "$OUTPUT_R2"

        echo "$SAMPLE,EXTRACTION_FAILED,0,0" >> "$SUMMARY"

        FAILED=$((FAILED + 1))

        continue

    fi

    ###########################################################################
    # COUNT EXTRACTED READS
    ###########################################################################

    if [[ -f "$OUTPUT_R1" ]]
    then
        R1_READS=$(awk 'END {print NR/4}' "$OUTPUT_R1")
    else
        R1_READS=0
    fi

    if [[ -f "$OUTPUT_R2" ]]
    then
        R2_READS=$(awk 'END {print NR/4}' "$OUTPUT_R2")
    else
        R2_READS=0
    fi

    ###########################################################################
    # HANDLE ZERO READS
    ###########################################################################

    if [[ "$R1_READS" -eq 0 ]]
    then

        echo
        echo "[WARNING] No reads found for:"
        echo "          $TAXON_NAME (taxid $TAXON_ID)"
        echo "          or any of its descendants."

        rm -f "$OUTPUT_R1" "$OUTPUT_R2"

        echo "$SAMPLE,NO_READS,0,0" >> "$SUMMARY"

        NO_READS=$((NO_READS + 1))

        continue

    fi

    ###########################################################################
    # SUCCESS
    ###########################################################################

    echo
    echo "[DONE] $SAMPLE"
    echo "R1 reads extracted : $R1_READS"
    echo "R2 reads extracted : $R2_READS"

    echo "$SAMPLE,DONE,$R1_READS,$R2_READS" >> "$SUMMARY"

    PROCESSED=$((PROCESSED + 1))

done

###############################################################################
# FINAL SUMMARY
###############################################################################

echo
echo "============================================================"
echo "             EXTRACTION COMPLETE"
echo "============================================================"

echo
echo "Target taxid  : $TAXON_ID"
echo "Target name   : $TAXON_NAME"

echo
echo "Input directory:"
echo "$INPUT_DIR"

echo
echo "Kraken directory:"
echo "$KRAKEN_DIR"

echo
echo "Output directory:"
echo "$OUTPUT_DIR"

echo
echo "Total R1 files found     : $FOUND"
echo "Successfully processed   : $PROCESSED"
echo "Already completed/skipped: $SKIPPED"
echo "No target reads          : $NO_READS"
echo "Missing R2               : $MISSING_R2"
echo "Failed                   : $FAILED"

echo
echo "Extracted reads:"
echo "$OUTPUT_DIR"

echo
echo "Summary:"
echo "$SUMMARY"

echo
echo "Finished: $(date)"

echo
echo "============================================================"