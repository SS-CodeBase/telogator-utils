#!/usr/bin/env bash

# Fail immediately if any command fails or if variables are unbound
set -euo pipefail

# ---- 1. Load configuration ----
# Use the provided config file, or default to 'template_4.tal.summary.config'
CONFIG="${1:-template_4.tal.summary.config}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found at: $CONFIG"
    exit 1
fi

source "$CONFIG"

# ---- 2. Validate required variables ----
# The script will abort here if the config file is missing any of these keys
: "${SPECIES_CHR_TABLE:?Missing SPECIES_CHR_TABLE in config}"
: "${RAW_DIR:?Missing RAW_DIR in config}"
: "${OUT_DIR:?Missing OUT_DIR in config}"
: "${SCRIPT:?Missing SCRIPT in config}"
: "${THREADS:?Missing THREADS in config}"
: "${ALLELE_MODE:=collapsed}" # Defaults to collapsed if not specified

export SCRIPT RAW_DIR OUT_DIR SPECIES_CHR_TABLE ALLELE_MODE

# ---- 3. Setup workspace ----
mkdir -p "$OUT_DIR/logs"

# Extract unique species/sample IDs from the 5th column of the reference table
echo "Extracting species list from reference table..."
mapfile -t SPECIES_LIST < <(cut -f5 "$SPECIES_CHR_TABLE" | tail -n +2 | sort -u)

echo "Found ${#SPECIES_LIST[@]} unique species/samples."
echo "Executing with $THREADS parallel threads (Allele mode: $ALLELE_MODE)..."

# ---- 4. Worker function for parallel execution ----
run_species() {
    local SPECIES="$1"
    local INPUT="${RAW_DIR}/${SPECIES}_tlens_by_allele.tsv"
    local OUTPUT="${OUT_DIR}/${SPECIES}.all.chr.tal.summary.${ALLELE_MODE}.tsv"
    local ERROR_FILE="${OUT_DIR}/${SPECIES}.error.txt"
    local LOG_FILE="${OUT_DIR}/logs/${SPECIES}.log"

    echo "[$SPECIES] Processing started" >> "$LOG_FILE"

    # Skip processing if input data is missing for this sample
    if [[ ! -f "$INPUT" ]]; then
        echo "$SPECIES: Input file missing at $INPUT" > "$ERROR_FILE"
        echo "[$SPECIES] Skipped - missing input" >> "$LOG_FILE"
        return 0
    fi

    # Execute Python script and redirect standard output and errors to the log
    if ! python3 "$SCRIPT" \
        --input "$INPUT" \
        --species_table "$SPECIES_CHR_TABLE" \
        --species "$SPECIES" \
        --output "$OUTPUT" \
        --error_file "$ERROR_FILE" \
        --allele_mode "$ALLELE_MODE" >> "$LOG_FILE" 2>&1; then
        echo "[$SPECIES] ERROR during python execution. Check log." >> "$LOG_FILE"
    fi

    echo "[$SPECIES] Processing complete" >> "$LOG_FILE"
}

# Export the function so xargs can see it in child subshells
export -f run_species

# ---- 5. Execute in parallel ----
# Passes the array of species into xargs, running $THREADS jobs simultaneously
printf "%s\n" "${SPECIES_LIST[@]}" | xargs -I {} -P "$THREADS" bash -c 'run_species "$@"' _ {}

echo "Pipeline execution completed for all samples."
