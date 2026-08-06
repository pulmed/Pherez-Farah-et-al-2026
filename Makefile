R_ENV ?= pherez-farah-r_env

.PHONY: validate smoke parse shell-check ensure-layers trajectory-smoke cellchat-smoke biotin-qc cluster7-heatmaps biotin-cellchat-panels pseudobulk

validate:
	conda run -n $(R_ENV) Rscript scripts/validate_repository_inputs.R

smoke:
	R_ENV=$(R_ENV) scripts/run_smoke_tests.sh

parse:
	conda run -n $(R_ENV) Rscript --vanilla -e 'source("renv/activate.R"); invisible(lapply(list.files("scripts", pattern="[.]R$$", recursive=TRUE, full.names=TRUE), parse))'

shell-check:
	bash -n scripts/build_object/cellranger_multi.sh scripts/qc/scan_contaminants_per_run.sh scripts/run_smoke_tests.sh

ensure-layers:
	conda run -n $(R_ENV) Rscript scripts/build_object/ensure_analysis_layers.R

trajectory-smoke:
	SMOKE_TEST=1 conda run -n $(R_ENV) Rscript scripts/trajectory/run_trajectory_analysis.R

cellchat-smoke:
	SMOKE_TEST=1 conda run -n $(R_ENV) Rscript scripts/cellchat/run_biotin_comparisons.R

biotin-qc:
	conda run -n $(R_ENV) Rscript scripts/qc/biotin_thresholds_and_status.R

cluster7-heatmaps:
	conda run -n $(R_ENV) Rscript scripts/visualization/cluster7_marker_heatmaps.R

biotin-cellchat-panels:
	conda run -n $(R_ENV) Rscript scripts/visualization/biotin_cellchat_feature_panels.R

pseudobulk:
	conda run -n $(R_ENV) Rscript scripts/pseudobulk/run_cluster_pseudobulk_dge.R
