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
export RENV_CONFIG_SYNCHRONIZED_CHECK="${RENV_CONFIG_SYNCHRONIZED_CHECK:-false}"

echo "[1/10] Parsing R scripts"
conda run -n "${R_ENV}" Rscript --vanilla -e 'invisible(lapply(list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE), parse))'

echo "[2/10] Checking shell syntax"
bash -n scripts/build_object/cellranger_multi.sh scripts/qc/scan_contaminants_per_run.sh scripts/run_smoke_tests.sh

echo "[3/10] Validating canonical Seurat object"
conda run -n "${R_ENV}" Rscript scripts/validate_repository_inputs.R

echo "[4/10] Human GSE221553 object builder smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/build_object/build_human_gse221553_seurat_object.R

echo "[5/10] Trajectory smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/trajectory/run_trajectory_analysis.R

echo "[6/10] CellChat smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/cellchat/run_biotin_comparisons.R

echo "[7/10] M1/M2 signature smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/signatures/compare_m1_m2_signatures_with_violins.R

echo "[8/10] Cross-presentation signature smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/signatures/score_cross_presentation_signature.R

echo "[9/10] Mouse-to-human ortholog helper smoke test"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/signatures/cross_species/translate_mouse_signatures_to_human.R

echo "[10/10] Human script smoke tests"
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/qc/human/global_delta_min5_cell_qc.R
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/signatures/human/score_monocyte_classical_nonclassical_signatures.R
SMOKE_TEST=1 conda run -n "${R_ENV}" Rscript scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R

echo "Smoke tests completed successfully."
