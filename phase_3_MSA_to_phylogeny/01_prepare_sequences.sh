#!/usr/bin/env bash

# ============================================================
# STEP 1: ORGANISE AND CLEAN SEQUENCES FOR MULTIPLE SEQUENCE
# ALIGNMENT
#
# Purpose:
#   1. Retain the original input FASTA unchanged
#   2. Remove duplicate FASTA headers
#   3. Preserve sequences with different headers, even when
#      their nucleotide sequences are identical
#   4. Convert nucleotide characters to uppercase
#   5. Preserve valid IUPAC nucleotide ambiguity codes
#   6. Replace invalid nucleotide characters with N
#   7. Normalise Windows/CRLF line endings to Unix format
#   8. Validate the final FASTA
#
# Output:
#   One final clean FASTA only.
#   No intermediate FASTA files are retained.
# ============================================================

set -euo pipefail

# ============================================================
# USER CONFIGURATION
# ============================================================

PROJECT_DIR="/home/fadele/ARS_NACA/ppr_phylo"

INPUT_FASTA="${PROJECT_DIR}/ppr_pub_study_og.fasta"

OUTPUT_DIR="${PROJECT_DIR}/01_sequence_preparation"

OUTPUT_FASTA="${OUTPUT_DIR}/ppr_pub_study_og_clean.fasta"

LOG_FILE="${OUTPUT_DIR}/sequence_cleaning.log"


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

mkdir -p "${OUTPUT_DIR}"


# ============================================================
# CHECK REQUIRED INPUT
# ============================================================

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Input FASTA not found:"
    echo "${INPUT_FASTA}"
    exit 1
fi

if [[ ! -s "${INPUT_FASTA}" ]]; then
    echo "ERROR: Input FASTA exists but is empty:"
    echo "${INPUT_FASTA}"
    exit 1
fi


# ============================================================
# CHECK REQUIRED COMMANDS
# ============================================================

for CMD in awk sed grep tr; do
    if ! command -v "${CMD}" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: ${CMD}"
        exit 1
    fi
done


# ============================================================
# START LOG
# ============================================================

{
    echo "============================================================"
    echo "PPR PHYLOGENY - STEP 1"
    echo "Sequence preparation and cleaning"
    echo "============================================================"
    echo
    echo "Date:        $(date)"
    echo "Project:     ${PROJECT_DIR}"
    echo "Input:       ${INPUT_FASTA}"
    echo "Output:      ${OUTPUT_FASTA}"
    echo
} > "${LOG_FILE}"


# ============================================================
# COUNT INPUT SEQUENCES
# ============================================================

INPUT_COUNT=$(grep -c '^>' "${INPUT_FASTA}" || true)

if [[ "${INPUT_COUNT}" -eq 0 ]]; then
    echo "ERROR: No FASTA headers were detected in the input file."
    exit 1
fi


# ============================================================
# CLEAN FASTA
#
# Rules:
#
# - FASTA headers are identified by "^>"
# - Header identity is used for deduplication
# - Only the FIRST occurrence of an identical header is retained
# - Different headers are always retained
# - Sequence characters are converted to uppercase
# - Valid IUPAC nucleotide characters are retained:
#
#   A C G T N R Y K M S W B D H V
#
# - Any other character is converted to N
# - Blank sequence lines are ignored
# - CR characters from Windows line endings are removed
#
# No intermediate FASTA is produced.
# ============================================================

awk '
BEGIN {
    FS = ""
    OFS = ""

    valid = "ACGTNR YKMSWBDHV"

    # Remove the accidental space from the definition.
    gsub(/ /, "", valid)

    current_header = ""
    skip = 0

    input_headers = 0
    duplicate_headers = 0
    invalid_characters = 0
    sequences_with_N = 0
    sequence_started = 0
}

# ------------------------------------------------------------
# FASTA HEADER
# ------------------------------------------------------------
/^>/ {

    # Remove carriage return from Windows/CRLF files.
    sub(/\r$/, "")

    input_headers++

    current_header = $0

    # Duplicate header:
    # skip the complete sequence belonging to this header.
    if (seen[current_header]++) {
        skip = 1
        duplicate_headers++
        next
    }

    skip = 0

    print current_header

    sequence_started = 0

    next
}

# ------------------------------------------------------------
# SEQUENCE
# ------------------------------------------------------------
{
    if (skip)
        next

    # Remove carriage return.
    sub(/\r$/, "")

    # Ignore completely blank lines.
    if ($0 == "")
        next

    sequence_started = 1

    line = toupper($0)
    cleaned = ""

    for (i = 1; i <= length(line); i++) {

        base = substr(line, i, 1)

        if (index(valid, base) > 0) {
            cleaned = cleaned base

            if (base == "N")
                has_N[current_header] = 1

        } else {
            cleaned = cleaned "N"
            invalid_characters++
        }
    }

    print cleaned
}

END {

    for (h in has_N)
        sequences_with_N++
}
' "${INPUT_FASTA}" > "${OUTPUT_FASTA}"


# ============================================================
# VERIFY OUTPUT WAS CREATED
# ============================================================

if [[ ! -s "${OUTPUT_FASTA}" ]]; then
    echo "ERROR: Clean FASTA was not generated successfully."
    rm -f "${OUTPUT_FASTA}"
    exit 1
fi


# ============================================================
# COUNT OUTPUT SEQUENCES
# ============================================================

OUTPUT_COUNT=$(grep -c '^>' "${OUTPUT_FASTA}" || true)


# ============================================================
# CHECK FOR DUPLICATE HEADERS IN OUTPUT
# ============================================================

DUPLICATE_HEADERS_REMAINING=$(
    grep '^>' "${OUTPUT_FASTA}" |
    sort |
    uniq -d |
    wc -l
)

if [[ "${DUPLICATE_HEADERS_REMAINING}" -ne 0 ]]; then
    echo "ERROR: Duplicate FASTA headers remain in the output."
    rm -f "${OUTPUT_FASTA}"
    exit 1
fi


# ============================================================
# CHECK FOR INVALID NUCLEOTIDE CHARACTERS
# ============================================================

INVALID_REMAINING=$(
    grep -v '^>' "${OUTPUT_FASTA}" |
    tr -d 'ACGTNR YKMSWBDHV\n\r' |
    wc -c
)

if [[ "${INVALID_REMAINING}" -ne 0 ]]; then
    echo "ERROR: Invalid nucleotide characters remain in output."
    rm -f "${OUTPUT_FASTA}"
    exit 1
fi


# ============================================================
# CHECK FOR WINDOWS CARRIAGE RETURNS
# ============================================================

CR_REMAINING=$(
    grep -Ul $'\r' "${OUTPUT_FASTA}" | wc -l
)

if [[ "${CR_REMAINING}" -ne 0 ]]; then
    echo "ERROR: Windows carriage returns remain in output."
    rm -f "${OUTPUT_FASTA}"
    exit 1
fi


# ============================================================
# CALCULATE DUPLICATES REMOVED
# ============================================================

DUPLICATES_REMOVED=$((INPUT_COUNT - OUTPUT_COUNT))


# ============================================================
# WRITE SUMMARY TO LOG
# ============================================================

{
    echo
    echo "-------------------- SUMMARY --------------------"
    echo
    echo "Input sequences:                  ${INPUT_COUNT}"
    echo "Output sequences:                 ${OUTPUT_COUNT}"
    echo "Duplicate headers removed:       ${DUPLICATES_REMOVED}"
    echo "Sequences containing N:           ${sequences_with_N:-0}"
    echo "Invalid characters replaced:      ${invalid_characters:-0}"
    echo "Duplicate headers remaining:      ${DUPLICATE_HEADERS_REMAINING}"
    echo "Invalid characters remaining:     ${INVALID_REMAINING}"
    echo "Windows CR characters remaining:  ${CR_REMAINING}"
    echo
    echo "Final FASTA:"
    echo "${OUTPUT_FASTA}"
    echo
    echo "Sequence preparation completed successfully."
    echo "--------------------------------------------------"
} >> "${LOG_FILE}"


# ============================================================
# DISPLAY SUMMARY
# ============================================================

cat "${LOG_FILE}"

echo
echo "============================================================"
echo "STEP 1 COMPLETED SUCCESSFULLY"
echo "============================================================"
echo
echo "Original input retained:"
echo "  ${INPUT_FASTA}"
echo
echo "Clean FASTA:"
echo "  ${OUTPUT_FASTA}"
echo
echo "Log:"
echo "  ${LOG_FILE}"
echo
