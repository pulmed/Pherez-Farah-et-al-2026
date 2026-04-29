#!/usr/bin/env bash
# ==============================================================================
# Script: 00_run_cellranger_multi.sh
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Run Cell Ranger multi for GEX, antibody-capture, and V(D)J libraries.
# Inputs: Raw FASTQ folders, 10x reference directories, and feature-reference CSV.
# Outputs: Cell Ranger multi output directories and generated config CSV files.
# Dependencies: bash, find, sed, Cell Ranger 9.0.1.
# Notes: Samples are auto-detected from *_GEX folders in BASE_DIR.
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------
# Override these from the command line if needed, for example:
# BASE_DIR=/path/to/fastqs RESULTS_DIR=/path/to/results ./scripts/00_run_cellranger_multi.sh

BASE_DIR="${BASE_DIR:-/path/to/raw_fastqs}"
GEX_REF="${GEX_REF:-/path/to/refdata-gex}"
VDJ_REF="${VDJ_REF:-/path/to/refdata-vdj}"
FEATURE_REF="${FEATURE_REF:-/path/to/feature_reference.csv}"
RESULTS_DIR="${RESULTS_DIR:-/path/to/results}"

mkdir -p "${RESULTS_DIR}"

# ------------------------------------------------------------------------------
# 2. FUNCTIONS
# ------------------------------------------------------------------------------

get_sample_id() {
    local fastq_dir="$1"

    find "${fastq_dir}" -name "*.fastq.gz" \
        | head -n 1 \
        | sed -E 's|.*/||' \
        | sed -E 's|_S[0-9]+_L.*||'
}

# ------------------------------------------------------------------------------
# 3. AUTO-DETECT SAMPLES
# ------------------------------------------------------------------------------

SAMPLES="$(find "${BASE_DIR}" -maxdepth 1 -type d -name "*_GEX" \
    | sed -E 's|.*/||' \
    | sed -E 's|_GEX$||' \
    | sort)"

if [[ -z "${SAMPLES}" ]]; then
    echo "[!] No samples detected in ${BASE_DIR}"
    echo "[!] Expected folders named: <sample>_GEX, <sample>_CSP, <sample>_VDJ"
    exit 1
fi

# ------------------------------------------------------------------------------
# 4. RUN CELL RANGER MULTI
# ------------------------------------------------------------------------------

for SAMPLE in ${SAMPLES}; do
    echo "[→] Processing ${SAMPLE}"

    GEX_FASTQ="${BASE_DIR}/${SAMPLE}_GEX"
    CSP_FASTQ="${BASE_DIR}/${SAMPLE}_CSP"
    VDJ_FASTQ="${BASE_DIR}/${SAMPLE}_VDJ"

    if [[ ! -d "${GEX_FASTQ}" ]]; then
        echo "[!] Missing GEX folder for ${SAMPLE}, skipping"
        continue
    fi

    if [[ ! -d "${CSP_FASTQ}" ]]; then
        echo "[!] Missing CSP folder for ${SAMPLE}, skipping"
        continue
    fi

    if [[ ! -d "${VDJ_FASTQ}" ]]; then
        echo "[!] Missing VDJ folder for ${SAMPLE}, skipping"
        continue
    fi

    GEX_ID="$(get_sample_id "${GEX_FASTQ}")"
    CSP_ID="$(get_sample_id "${CSP_FASTQ}")"
    VDJ_ID="$(get_sample_id "${VDJ_FASTQ}")"

    if [[ -z "${GEX_ID}" || -z "${CSP_ID}" || -z "${VDJ_ID}" ]]; then
        echo "[!] Could not detect FASTQ IDs for ${SAMPLE}, skipping"
        continue
    fi

    CONFIG_FILE="${RESULTS_DIR}/${SAMPLE}_multi_config.csv"

    cat > "${CONFIG_FILE}" <<EOF
[gene-expression]
reference,${GEX_REF}
create-bam,true

[vdj]
reference,${VDJ_REF}

[libraries]
fastq_id,fastqs,feature_types
${GEX_ID},${GEX_FASTQ},Gene Expression
${CSP_ID},${CSP_FASTQ},Antibody Capture
${VDJ_ID},${VDJ_FASTQ},VDJ-T

[feature]
reference,${FEATURE_REF}
EOF

    echo "[→] Running cellranger multi for ${SAMPLE}"

    cellranger multi \
        --id="${SAMPLE}_multi" \
        --csv="${CONFIG_FILE}" \
        --output-dir="${RESULTS_DIR}/${SAMPLE}_multi"

    echo "[✓] ${SAMPLE} done"
done