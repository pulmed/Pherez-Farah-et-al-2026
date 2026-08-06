#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: scripts/run_smoke_tests.sh
# Original file: repository smoke-test utility
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Run lightweight repository checks without executing expensive full analyses.
# Inputs: Repository scripts and canonical Seurat object, default data/seurat.rds.
# Outputs: Console smoke-test report.
# Dependencies: bash, Rscript, conda environment pherez-farah-r_env.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Full CellChat and trajectory analyses are skipped through SMOKE_TEST=1.
# ------------------------------------------------------------------------------

set -euo pipefail

R_ENV="${R_ENV:-pherez-farah-r_env}"

echo "[1/5] Parsing R scripts"
conda run -n "${R_ENV}" Rscript --vanilla -e 'source("renv/activate.R"); invisible(lapply(list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE), parse))'

echo "[2/5] Checking shell syntax"
bash -n scripts/build_object/cellranger_multi.sh scripts/qc/scan_contaminants_per_run.sh scripts/run_smoke_tests.sh

echo "[3/5] Validating canonical Seurat object"
conda run -n "${R_ENV}" Rscript scripts/validate_repository_inputs.R

echo "[4/5] Trajectory smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/trajectory/run_trajectory_analysis.R

echo "[5/5] CellChat smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/cellchat/run_biotin_comparisons.R

echo "Smoke tests completed successfully."
