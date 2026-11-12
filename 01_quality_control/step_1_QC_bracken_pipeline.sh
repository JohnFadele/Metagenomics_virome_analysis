#!/usr/bin/env bash

# ==============================================================================
# UNIFIED VIRAL METAGENOMICS PIPELINE
# Steps: Rename -> fastp QC -> Sequential Host Depletion -> Kraken2 -> Bracken
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# ==============================================================================
# 1. USER CONFIGURATION
# Configurable paths, thread counts, parameters, and stage toggles.
# ==============================================================================

# Directory Setup
INPUT_DIR="/home/fadele/ppr_extra"                      # Raw paired FASTQ files
OUTPUT_BASE="/home/fadele/ppr_extra/results"           # Parent output directory
HOST_DB_DIR="/home/fadele/host_genomes/domestic"  # Directory containing host Bowtie2 indexes
KRAKEN_DB="/home/fadele/Kraken_full/k2_standard_20260626" # Kraken2 database

# Stage Enable/Disable Switches (true/false)
ENABLE_RENAME_CLEANUP=true   # Strips Illumina '_001' suffixes
ENABLE_QC_TRIMMING=true      # Runs fastp trimming & QC
ENABLE_HOST_DEPLETION=true   # Sequentially maps/filters host reads
ENABLE_KRAKEN2=true          # Taxonomic profiling
ENABLE_BRACKEN=true          # Abundance re-estimation

# General Performance
THREADS=8

# fastp Settings
CUT_WINDOW_SIZE=4
CUT_MEAN_QUALITY=20
MIN_READ_LENGTH=50

# Kraken2 Settings
KRAKEN_CONFIDENCE=0.0        # Precision threshold (0.0 recommended)
KRAKEN_MIN_HIT_GROUPS=3      # K-mer cluster threshold to reduce false positives

# Bracken Settings
BRACKEN_READ_LEN=150         # Read length for model estimation
BRACKEN_LEVEL="S"            # Taxonomic level (S=Species, G=Genus, etc.)

# ==============================================================================
# 2. AUTOMATIC DIRECTORY ROUTING
# Derived directories created dynamically based on OUTPUT_BASE.
# ==============================================================================

TRIM_DIR="${OUTPUT_BASE}/fastp_trimmed"
REPORT_DIR="${OUTPUT_BASE}/fastp_reports"
LOG_DIR="${OUTPUT_BASE}/logs"
DEHOST_DIR="${OUTPUT_BASE}/dehosted"
KRAKEN_DIR="${OUTPUT_BASE}/kraken2"
BRACKEN_DIR="${OUTPUT_BASE}/bracken"

QC_SUMMARY_CSV="${OUTPUT_BASE}/fastp_QC_summary.csv"
DEHOST_SUMMARY="${DEHOST_DIR}/dehost_summary.csv"
DEHOST_FAILED="${DEHOST_DIR}/failed_samples.txt"
ERROR_LOG="${LOG_DIR}/pipeline_errors.log"

mkdir -p "$OUTPUT_BASE" "$TRIM_DIR" "$REPORT_DIR" "$LOG_DIR" \
         "$DEHOST_DIR" "$KRAKEN_DIR" "$BRACKEN_DIR"

echo "=============================================================================="
echo "STARTING METAGENOMICS PIPELINE"
echo "Date: $(date)"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_BASE"
echo "=============================================================================="

# Check system dependencies
for tool in fastp jq bowtie2 kraken2 bracken; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[ERROR] Required tool '$tool' is missing from PATH."
        exit 1
    fi
done

# ==============================================================================
# STEP 0: RENAME FASTQ FILES (Remove _001 Suffix)
# ==============================================================================
if [[ "$ENABLE_RENAME_CLEANUP" == true ]]; then
    echo
    echo "--- [STEP 0] Renaming FASTQ files (Stripping _001 suffix) ---"
    for f in "$INPUT_DIR"/*_R[12]_001.fastq.gz "$INPUT_DIR"/*_R[12]_001.fastq; do
        [[ -f "$f" ]] || continue
        new="${f/_001./.}"
        echo "  Renaming: $(basename "$f") -> $(basename "$new")"
        mv -- "$f" "$new"
    done
fi

# ==============================================================================
# STEP 1: FASTP QUALITY CONTROL & TRIMMING
# ==============================================================================
if [[ "$ENABLE_QC_TRIMMING" == true ]]; then
    echo
    echo "--- [STEP 1] Running fastp Quality Control and Adapter Trimming ---"

    for R1 in "$INPUT_DIR"/*_R1.fastq "$INPUT_DIR"/*_R1.fastq.gz; do
        [[ -f "$R1" ]] || continue

        if [[ "$R1" == *.fastq.gz ]]; then
            R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        else
            R2="${R1/_R1.fastq/_R2.fastq}"
        fi

        SAMPLE=$(basename "$R1")
        SAMPLE="${SAMPLE%_R1.fastq.gz}"
        SAMPLE="${SAMPLE%_R1.fastq}"

        if [[ ! -f "$R2" ]]; then
            echo "[WARNING] Missing R2 for $SAMPLE. Skipping."
            continue
        fi

        TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
        TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
        HTML_REP="${REPORT_DIR}/${SAMPLE}.fastp.html"
        JSON_REP="${REPORT_DIR}/${SAMPLE}.fastp.json"

        # Resume support
        if [[ -s "$TRIM_R1" && -s "$TRIM_R2" && -s "$HTML_REP" && -s "$JSON_REP" ]]; then
            echo "[SKIP] $SAMPLE fastp QC already completed."
            continue
        fi

        echo "[QC] Processing: $SAMPLE"
        fastp \
            --in1 "$R1" --in2 "$R2" \
            --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
            --detect_adapter_for_pe \
            --cut_front --cut_tail \
            --cut_window_size "$CUT_WINDOW_SIZE" \
            --cut_mean_quality "$CUT_MEAN_QUALITY" \
            --length_required "$MIN_READ_LENGTH" \
            --thread "$THREADS" \
            --html "$HTML_REP" --json "$JSON_REP"
    done

    # Generate Aggregated QC Report
    echo "Sample,Reads_before,Reads_after,Reads_retained(%),Bases_before,Bases_after,Q20_before(%),Q30_before(%),Q20_after(%),Q30_after(%),GC_before(%),GC_after(%),Duplication_rate(%)" > "$QC_SUMMARY_CSV"
    
    for JSON in "$REPORT_DIR"/*.fastp.json; do
        [[ -s "$JSON" ]] || continue
        SNAME=$(basename "$JSON" .fastp.json)
        
        jq -r '
        [
            .summary.before_filtering.total_reads,
            .summary.after_filtering.total_reads,
            (if .summary.before_filtering.total_reads > 0 then (.summary.after_filtering.total_reads / .summary.before_filtering.total_reads * 100) else 0 end),
            .summary.before_filtering.total_bases,
            .summary.after_filtering.total_bases,
            (.summary.before_filtering.q20_rate * 100),
            (.summary.before_filtering.q30_rate * 100),
            (.summary.after_filtering.q20_rate * 100),
            (.summary.after_filtering.q30_rate * 100),
            (.summary.before_filtering.gc_content * 100),
            (.summary.after_filtering.gc_content * 100),
            (.duplication.rate * 100)
        ] | @csv' "$JSON" | sed "s/^/\"$SNAME\",/" >> "$QC_SUMMARY_CSV"
    done
    echo "[DONE] Aggregated QC report saved to: $QC_SUMMARY_CSV"
fi

# ==============================================================================
# STEP 2: SEQUENTIAL MULTI-HOST DEPLETION
# ==============================================================================
if [[ "$ENABLE_HOST_DEPLETION" == true ]]; then
    echo
    echo "--- [STEP 2] Running Sequential Host Depletion ---"

    echo "Sample,Status" > "$DEHOST_SUMMARY"
    > "$DEHOST_FAILED"

    # Detect indexed host genomes
    HOSTS=()
    for FASTA in "$HOST_DB_DIR"/*.fa "$HOST_DB_DIR"/*.fasta; do
        [[ -f "$FASTA" ]] || continue
        BASE="${FASTA%.*}"
        if [[ -f "${BASE}.1.bt2" || -f "${BASE}.1.bt2l" ]]; then
            HOSTS+=("$(basename "$BASE")")
        fi
    done

    if [[ ${#HOSTS[@]} -eq 0 ]]; then
        echo "[WARNING] No Bowtie2 index files detected in $HOST_DB_DIR. Skipping host depletion."
    else
        echo "Detected ${#HOSTS[@]} host indexes: ${HOSTS[*]}"

        for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
            [[ -f "$R1" ]] || continue
            SAMPLE=$(basename "$R1" "_R1.trimmed.fastq.gz")
            R2="$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz"

            FINAL_R1="$DEHOST_DIR/${SAMPLE}_R1.clean.fastq.gz"
            FINAL_R2="$DEHOST_DIR/${SAMPLE}_R2.clean.fastq.gz"

            if [[ -f "$FINAL_R1" && -f "$FINAL_R2" ]]; then
                echo "[SKIP] $SAMPLE dehosting already completed."
                echo "$SAMPLE,SKIPPED" >> "$DEHOST_SUMMARY"
                continue
            fi

            CURRENT_R1="$R1"
            CURRENT_R2="$R2"
            SUCCESS=1

            for HOST in "${HOSTS[@]}"; do
                echo "[HOST DEPLETION] Removing $HOST from $SAMPLE"
                PREFIX="$DEHOST_DIR/${SAMPLE}_${HOST}"

                if bowtie2 \
                    --sensitive \
                    --end-to-end \
                    -p "$THREADS" \
                    -x "$HOST_DB_DIR/$HOST" \
                    -1 "$CURRENT_R1" \
                    -2 "$CURRENT_R2" \
                    --un-conc-gz "${PREFIX}_%.fastq.gz" \
                    -S /dev/null \
                    2>> "$LOG_DIR/${SAMPLE}_dehost.log"; then

                    NEXT_R1="${PREFIX}_1.fastq.gz"
                    NEXT_R2="${PREFIX}_2.fastq.gz"

                    if [[ "$CURRENT_R1" != "$R1" ]]; then
                        rm -f "$CURRENT_R1" "$CURRENT_R2"
                    fi
                    CURRENT_R1="$NEXT_R1"
                    CURRENT_R2="$NEXT_R2"
                else
                    echo "[ERROR] Host depletion failed for $SAMPLE on $HOST"
                    SUCCESS=0
                    break
                fi
            done

            if [[ $SUCCESS -eq 1 && -f "$CURRENT_R1" && -f "$CURRENT_R2" ]]; then
                mv "$CURRENT_R1" "$FINAL_R1"
                mv "$CURRENT_R2" "$FINAL_R2"
                echo "$SAMPLE,DONE" >> "$DEHOST_SUMMARY"
            else
                echo "$SAMPLE,FAILED" >> "$DEHOST_SUMMARY"
                echo "$SAMPLE" >> "$DEHOST_FAILED"
                rm -f "$DEHOST_DIR/${SAMPLE}_"*_1.fastq.gz "$DEHOST_DIR/${SAMPLE}_"*_2.fastq.gz
            fi
        done
    fi
fi

# ==============================================================================
# STEP 3: KRAKEN2 TAXONOMIC CLASSIFICATION
# ==============================================================================
if [[ "$ENABLE_KRAKEN2" == true ]]; then
    echo
    echo "--- [STEP 3] Running Kraken2 Taxonomic Classification ---"

    # Select dehosted directory if enabled and non-empty, otherwise use trimmed reads
    if [[ "$ENABLE_HOST_DEPLETION" == true && -n "$(ls -A "$DEHOST_DIR"/*_R1.clean.fastq.gz 2>/dev/null)" ]]; then
        KRAKEN_INPUT_DIR="$DEHOST_DIR"
        SUFFIX_R1="_R1.clean.fastq.gz"
        SUFFIX_R2="_R2.clean.fastq.gz"
    else
        KRAKEN_INPUT_DIR="$TRIM_DIR"
        SUFFIX_R1="_R1.trimmed.fastq.gz"
        SUFFIX_R2="_R2.trimmed.fastq.gz"
    fi

    for R1 in "$KRAKEN_INPUT_DIR"/*"$SUFFIX_R1"; do
        [[ -f "$R1" ]] || continue
        SAMPLE=$(basename "$R1" "$SUFFIX_R1")
        R2="$KRAKEN_INPUT_DIR/${SAMPLE}${SUFFIX_R2}"

        OUTPUT="$KRAKEN_DIR/${SAMPLE}.kraken.output"
        REPORT="$KRAKEN_DIR/${SAMPLE}.report.txt"

        if [[ -s "$OUTPUT" && -s "$REPORT" ]]; then
            echo "[SKIP] $SAMPLE Kraken2 profiling already completed."
            continue
        fi

        echo "[KRAKEN2] Processing: $SAMPLE"
        kraken2 \
            --db "$KRAKEN_DB" \
            --threads "$THREADS" \
            --paired \
            --gzip-compressed \
            --confidence "$KRAKEN_CONFIDENCE" \
            --minimum-hit-groups "$KRAKEN_MIN_HIT_GROUPS" \
            --output "$OUTPUT" \
            --report "$REPORT" \
            "$R1" "$R2"
    done
fi

# ==============================================================================
# STEP 4: BRACKEN SPECIES-LEVEL RE-ESTIMATION
# ==============================================================================
if [[ "$ENABLE_BRACKEN" == true ]]; then
    echo
    echo "--- [STEP 4] Running Bracken Abundance Re-estimation ---"

    for REPORT in "$KRAKEN_DIR"/*.report.txt; do
        [[ -f "$REPORT" ]] || continue
        SAMPLE=$(basename "$REPORT" .report.txt)

        OUTFILE="$BRACKEN_DIR/${SAMPLE}.species.bracken.txt"
        REPORTFILE="$BRACKEN_DIR/${SAMPLE}.species.report.txt"

        if [[ -f "$OUTFILE" && -f "$REPORTFILE" ]]; then
            echo "[SKIP] $SAMPLE Bracken re-estimation already completed."
            continue
        fi

        echo "[BRACKEN] Processing: $SAMPLE"
        bracken \
            -d "$KRAKEN_DB" \
            -i "$REPORT" \
            -o "$OUTFILE" \
            -w "$REPORTFILE" \
            -r "$BRACKEN_READ_LEN" \
            -l "$BRACKEN_LEVEL" \
            -t "$THREADS"
    done
fi

echo
echo "=============================================================================="
echo "PIPELINE EXECUTED SUCCESSFULLY"
echo "Finished at: $(date)"
echo "Output Directory Structure:"
echo "  QC & Reports:       $REPORT_DIR"
echo "  De-hosted Reads:    $DEHOST_DIR"
echo "  Kraken2 Profiles:   $KRAKEN_DIR"
echo "  Bracken Estimates:  $BRACKEN_DIR"
echo "=============================================================================="