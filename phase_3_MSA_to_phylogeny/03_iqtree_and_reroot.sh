#!/bin/bash

# ============================================================
# STEP 4: IQ-TREE PHYLOGENY + TREE REROOTING
#
# Workflow:
#   1. Locate *_trimmed.fasta from Step 2
#   2. Build maximum-likelihood tree with IQ-TREE
#   3. Generate .treefile and .contree
#   4. Reroot both trees using the specified outgroup
#   5. Retain both original and rooted trees
#
# ============================================================

set -Eeuo pipefail

# ============================================================
# CONFIGURATION
# ============================================================

# Change this for a new project
PROJECT_DIR="/home/fadele/ARS_NACA/ppr_phylo"

# Outgroup accession for this analysis
OUTGROUP="NC_006296.2"

# Number of CPU threads
THREADS=8

# Input directory containing the trimmed FASTA
ALIGNMENT_DIR="${PROJECT_DIR}/02_alignment"

# Output directory
OUTPUT_DIR="${PROJECT_DIR}/03_IQTREE"

# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

mkdir -p "${OUTPUT_DIR}"

# ============================================================
# CHECK REQUIRED PROGRAMS
# ============================================================

for PROGRAM in iqtree nw_reroot; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} was not found in PATH."
        echo "Please activate the appropriate conda environment."
        exit 1
    fi
done

# ============================================================
# LOCATE TRIMMED FASTA
# ============================================================

mapfile -t TRIMMED_FILES < <(
    find "${ALIGNMENT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*_trimmed.fasta" \
        -print | sort
)

if [[ ${#TRIMMED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No *_trimmed.fasta file was found in:"
    echo "${ALIGNMENT_DIR}"
    exit 1
fi

if [[ ${#TRIMMED_FILES[@]} -gt 1 ]]; then
    echo "ERROR: More than one *_trimmed.fasta file was found:"
    printf '%s\n' "${TRIMMED_FILES[@]}"
    echo
    echo "The script expects one trimmed FASTA in:"
    echo "${ALIGNMENT_DIR}"
    exit 1
fi

INPUT_FILE="${TRIMMED_FILES[0]}"

# ============================================================
# DETERMINE ANALYSIS NAME
# ============================================================

INPUT_BASENAME="$(basename "${INPUT_FILE}")"

# Remove .fasta
ANALYSIS_NAME="${INPUT_BASENAME%.fasta}"

# Remove _trimmed
ANALYSIS_NAME="${ANALYSIS_NAME%_trimmed}"

# ============================================================
# DEFINE TREE PREFIX
# ============================================================

TREE_PREFIX="${OUTPUT_DIR}/${ANALYSIS_NAME}"

# Original IQ-TREE outputs
TREEFILE="${TREE_PREFIX}.treefile"
CONTREE="${TREE_PREFIX}.contree"

# Rooted outputs
ROOTED_TREEFILE="${OUTPUT_DIR}/${ANALYSIS_NAME}_rooted.treefile"
ROOTED_CONTREE="${OUTPUT_DIR}/${ANALYSIS_NAME}_rooted.contree"

# Log
LOG_FILE="${OUTPUT_DIR}/${ANALYSIS_NAME}_phylogeny.log"

# ============================================================
# START LOG
# ============================================================

{
    echo "============================================================"
    echo "IQ-TREE PHYLOGENY AND REROOTING"
    echo "============================================================"
    echo "Start time:       $(date)"
    echo
    echo "Project directory:"
    echo "  ${PROJECT_DIR}"
    echo
    echo "Input FASTA:"
    echo "  ${INPUT_FILE}"
    echo
    echo "Analysis name:"
    echo "  ${ANALYSIS_NAME}"
    echo
    echo "Threads:"
    echo "  ${THREADS}"
    echo
    echo "Outgroup:"
    echo "  ${OUTGROUP}"
    echo
    echo "Output directory:"
    echo "  ${OUTPUT_DIR}"
    echo "============================================================"
    echo
} | tee "${LOG_FILE}"

# ============================================================
# STEP 1 — BUILD MAXIMUM-LIKELIHOOD TREE
# ============================================================

echo "STEP 1: Building maximum-likelihood tree..." | tee -a "${LOG_FILE}"
echo | tee -a "${LOG_FILE}"

iqtree \
    -s "${INPUT_FILE}" \
    -pre "${TREE_PREFIX}" \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt "${THREADS}" \
    2>&1 | tee -a "${LOG_FILE}"

# ============================================================
# CHECK IQ-TREE OUTPUTS
# ============================================================

if [[ ! -s "${TREEFILE}" ]]; then
    echo "ERROR: IQ-TREE did not produce:"
    echo "${TREEFILE}"
    exit 1
fi

if [[ ! -s "${CONTREE}" ]]; then
    echo "ERROR: IQ-TREE did not produce:"
    echo "${CONTREE}"
    exit 1
fi

echo | tee -a "${LOG_FILE}"
echo "IQ-TREE successfully produced:" | tee -a "${LOG_FILE}"
echo "  ${TREEFILE}" | tee -a "${LOG_FILE}"
echo "  ${CONTREE}" | tee -a "${LOG_FILE}"

# ============================================================
# STEP 2 — VERIFY OUTGROUP IN BOTH TREES
# ============================================================

echo | tee -a "${LOG_FILE}"
echo "STEP 2: Checking outgroup..." | tee -a "${LOG_FILE}"
echo "Outgroup: ${OUTGROUP}" | tee -a "${LOG_FILE}"

if ! grep -Fq "${OUTGROUP}" "${TREEFILE}"; then
    echo "ERROR: Outgroup '${OUTGROUP}' was not found in:"
    echo "${TREEFILE}"
    exit 1
fi

if ! grep -Fq "${OUTGROUP}" "${CONTREE}"; then
    echo "ERROR: Outgroup '${OUTGROUP}' was not found in:"
    echo "${CONTREE}"
    exit 1
fi

echo "Outgroup found in both trees." | tee -a "${LOG_FILE}"

# ============================================================
# STEP 3 — REROOT THE ML TREE
# ============================================================

echo | tee -a "${LOG_FILE}"
echo "STEP 3: Rerooting maximum-likelihood tree..." | tee -a "${LOG_FILE}"

nw_reroot \
    "${TREEFILE}" \
    "${OUTGROUP}" \
    > "${ROOTED_TREEFILE}"

if [[ ! -s "${ROOTED_TREEFILE}" ]]; then
    echo "ERROR: Rooted .treefile was not created."
    exit 1
fi

echo "Rooted ML tree created:" | tee -a "${LOG_FILE}"
echo "  ${ROOTED_TREEFILE}" | tee -a "${LOG_FILE}"

# ============================================================
# STEP 4 — REROOT THE CONSENSUS TREE
# ============================================================

echo | tee -a "${LOG_FILE}"
echo "STEP 4: Rerooting consensus tree..." | tee -a "${LOG_FILE}"

nw_reroot \
    "${CONTREE}" \
    "${OUTGROUP}" \
    > "${ROOTED_CONTREE}"

if [[ ! -s "${ROOTED_CONTREE}" ]]; then
    echo "ERROR: Rooted .contree was not created."
    exit 1
fi

echo "Rooted consensus tree created:" | tee -a "${LOG_FILE}"
echo "  ${ROOTED_CONTREE}" | tee -a "${LOG_FILE}"

# ============================================================
# FINAL SUMMARY
# ============================================================

{
    echo
    echo "============================================================"
    echo "PHYLOGENETIC ANALYSIS COMPLETED SUCCESSFULLY"
    echo "============================================================"
    echo "End time: $(date)"
    echo
    echo "INPUT"
    echo "  ${INPUT_FILE}"
    echo
    echo "OUTGROUP"
    echo "  ${OUTGROUP}"
    echo
    echo "ORIGINAL UNROOTED ML TREE"
    echo "  ${TREEFILE}"
    echo
    echo "ORIGINAL UNROOTED CONSENSUS TREE"
    echo "  ${CONTREE}"
    echo
    echo "ROOTED ML TREE"
    echo "  ${ROOTED_TREEFILE}"
    echo
    echo "ROOTED CONSENSUS TREE"
    echo "  ${ROOTED_CONTREE}"
    echo
    echo "============================================================"
    echo "The rooted .treefile is the primary tree for visualisation."
    echo "============================================================"

} | tee -a "${LOG_FILE}"

echo
echo "SUCCESS: IQ-TREE construction and tree rerooting are complete."
