#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: scripts/qc/scan_contaminants_per_run.sh
# Original file: contamination scan shell workflow
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Identify reads aligning to a contaminant FASTA and count contaminated cell barcodes per run.
# Inputs: Contaminant FASTA and Cell Ranger sample_alignments.bam files.
# Outputs: Per-run FASTQ files, contaminant SAM hits, matched read names, contaminant read SAM lines, and contaminated cell-barcode count tables.
# Dependencies: bash, bowtie2, samtools, awk, sort, grep, cut.
# Environment: contam_scan conda environment.
# Notes:
# - Run with the contam_scan conda environment from envs/contam_scan_env.yml.
# ------------------------------------------------------------------------------

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------
CONTAM_FASTA="${CONTAM_FASTA:-data/contaminants.fa}"
CELLRANGER_RESULTS_DIR="${CELLRANGER_RESULTS_DIR:-output/cellranger_multi}"
OUTDIR="${OUTDIR:-output/contam_per_run}"
RUNS=(R1 R2 R3 R4 R5 R6)
THREADS=4

# ------------------------------------------------------------------------------
# 2. CHECK INPUTS AND CREATE OUTPUT DIRECTORY
# ------------------------------------------------------------------------------
if [[ ! -f "${CONTAM_FASTA}" ]]; then
  echo "ERROR: Contaminant FASTA not found: ${CONTAM_FASTA}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

# ------------------------------------------------------------------------------
# 3. BUILD CONTAMINANT INDEX
# ------------------------------------------------------------------------------
INDEX_PREFIX="${OUTDIR}/contaminant_index"

echo "[->] Building contaminant bowtie2 index"
bowtie2-build "${CONTAM_FASTA}" "${INDEX_PREFIX}" >/dev/null

# ------------------------------------------------------------------------------
# 4. PROCESS EACH RUN
# ------------------------------------------------------------------------------
for RUN in "${RUNS[@]}"; do
  echo "[->] Processing ${RUN}"

  BAM="${CELLRANGER_RESULTS_DIR}/${RUN}_multi/outs/per_sample_outs/${RUN}_multi/count/sample_alignments.bam"
  RUN_DIR="${OUTDIR}/${RUN}"

  if [[ ! -f "${BAM}" ]]; then
    echo "WARNING: BAM not found for ${RUN}: ${BAM}. Skipping." >&2
    continue
  fi

  mkdir -p "${RUN_DIR}"

  # Extract paired reads from Cell Ranger BAM.
  samtools fastq \
    -@ "${THREADS}" \
    -1 "${RUN_DIR}/R1.fastq" \
    -2 "${RUN_DIR}/R2.fastq" \
    -0 /dev/null \
    -s /dev/null \
    -n \
    "${BAM}"

  # Align extracted reads to the contaminant reference.
  bowtie2 \
    -x "${INDEX_PREFIX}" \
    -1 "${RUN_DIR}/R1.fastq" \
    -2 "${RUN_DIR}/R2.fastq" \
    -S "${RUN_DIR}/hits.sam" \
    -p "${THREADS}" \
    --quiet

  # Collect unique read names that aligned to the contaminant reference.
  awk '$1 !~ /^@/ {print $1}' "${RUN_DIR}/hits.sam" | sort -u > "${RUN_DIR}/matched_reads.txt"

  # Extract matching full SAM records from the original BAM.
  if [[ -s "${RUN_DIR}/matched_reads.txt" ]]; then
    samtools view "${BAM}" | grep -Ff "${RUN_DIR}/matched_reads.txt" > "${RUN_DIR}/contam_reads.sam" || true
  else
    : > "${RUN_DIR}/contam_reads.sam"
  fi

  # Extract and count Cell Ranger cell-barcode tags.
  grep -o 'CB:Z:[ACGT]*' "${RUN_DIR}/contam_reads.sam" \
    | cut -d: -f3 \
    | sort \
    | uniq -c \
    | sort -nr \
    > "${RUN_DIR}/contaminated_cells.txt" || true

  echo "[[OK]] ${RUN} done -> ${RUN_DIR}/contaminated_cells.txt"
done

echo "[[OK]] Contamination scan complete. Outputs written to ${OUTDIR}"
