#!/usr/bin/env bash

# ==============================================================================
# Autosomal TAL Aggregation Submitter
# ------------------------------------------------------------------------------
# This script reads a configuration file and executes the Python aggregation 
# script. It takes the individual summary TSV files generated in Step 1 and 
# compiles them into a single, multi-condition master matrix for the cohort.
# ==============================================================================

# Fail immediately if any command fails or if variables are unbound
set -euo pipefail

# ---- 1. Load configuration ----
# Use the provided config file, or default to 'template_4.tal.aggregate.config'
CONFIG_FILE="${1:-template_4.tal.aggregate.config}"

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Config file not found at: ${CONFIG_FILE}"
    exit 1
fi

# Load variables from the config file
# shellcheck disable=SC1090
source "${CONFIG_FILE}"

# ---- 2. Validate required variables ----
# The script will abort here if the config file is missing any of these keys.
# This prevents the Python script from failing halfway through due to missing paths.
: "${INPUT_DIR:?Missing INPUT_DIR in config}"
: "${OUTPUT_DIR:?Missing OUTPUT_DIR in config}"
: "${SCRIPT:?Missing SCRIPT in config}"
: "${FINAL_OUTPUT:?Missing FINAL_OUTPUT in config}"

echo "Starting TAL aggregation process..."
echo "Input directory: $INPUT_DIR"
echo "Output target: $FINAL_OUTPUT"

# ---- 3. Setup workspace ----
# Ensure the destination directory exists before writing to it
mkdir -p "$OUTPUT_DIR"

# ---- 4. Execute Aggregation ----
# Run the Python script using the validated paths
if python3 "$SCRIPT" \
    --input_dir "$INPUT_DIR" \
    --output "$FINAL_OUTPUT"; then
    
    echo "========================================"
    echo "Aggregation completed successfully!"
    echo "Results saved to: $FINAL_OUTPUT"
    echo "========================================"
else
    echo "ERROR: Aggregation failed during Python execution."
    exit 1
fi

exit 0
