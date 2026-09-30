R_ENV ?= pherez-farah-r_env
export RENV_CONFIG_SYNCHRONIZED_CHECK ?= false

.PHONY: validate smoke parse shell-check ensure-layers human-gse221553 human-rpca human-manuscript-object trajectory-smoke cellchat-smoke biotin-qc biotin-delta-sensitivity cluster7-heatmaps biotin-cellchat-panels pseudobulk-volcano-panels m1-m2 cross-presentation ucell-confirmation mouse-human-orthologues human-global-delta human-singler human-monocytes human-cluster2-isg pseudobulk

validate:
	conda run -n $(R_ENV) Rscript scripts/validate_repository_inputs.R

smoke:
	R_ENV=$(R_ENV) scripts/run_smoke_tests.sh

parse:
	conda run -n $(R_ENV) Rscript --vanilla -e 'invisible(lapply(list.files("scripts", pattern="[.]R$$", recursive=TRUE, full.names=TRUE), parse))'

shell-check:
	bash -n scripts/build_object/cellranger_multi.sh scripts/qc/scan_contaminants_per_run.sh scripts/run_smoke_tests.sh

ensure-layers:
	conda run -n $(R_ENV) Rscript scripts/build_object/ensure_analysis_layers.R

human-gse221553:
	conda run -n $(R_ENV) Rscript scripts/build_object/build_human_gse221553_seurat_object.R

human-rpca:
	conda run -n $(R_ENV) Rscript scripts/build_object/integrate_human_rpca_by_patient.R

human-manuscript-object:
	$(MAKE) human-gse221553
	$(MAKE) human-rpca

trajectory-smoke:
	SMOKE_TEST=1 conda run -n $(R_ENV) Rscript scripts/trajectory/run_trajectory_analysis.R

cellchat-smoke:
	SMOKE_TEST=1 conda run -n $(R_ENV) Rscript scripts/cellchat/run_biotin_comparisons.R

biotin-qc:
	conda run -n $(R_ENV) Rscript scripts/qc/biotin_thresholds_and_status.R

biotin-delta-sensitivity:
	conda run -n $(R_ENV) Rscript scripts/qc/biotin_delta_threshold_sensitivity.R

cluster7-heatmaps:
	conda run -n $(R_ENV) Rscript scripts/visualization/cluster7_marker_heatmaps.R

pseudobulk-volcano-panels:
	conda run -n $(R_ENV) Rscript scripts/visualization/plot_pseudobulk_volcano_panels.R

biotin-cellchat-panels:
	conda run -n $(R_ENV) Rscript scripts/visualization/biotin_cellchat_feature_panels.R

m1-m2:
	conda run -n $(R_ENV) Rscript scripts/signatures/compare_m1_m2_signatures_with_violins.R

cross-presentation:
	conda run -n $(R_ENV) Rscript scripts/signatures/score_cross_presentation_signature.R

ucell-confirmation:
	conda run -n $(R_ENV) Rscript scripts/signatures/confirm_ucell_signature_scores.R

mouse-human-orthologues:
	conda run -n $(R_ENV) Rscript scripts/signatures/cross_species/translate_mouse_signatures_to_human.R

human-global-delta:
	conda run -n $(R_ENV) Rscript scripts/qc/human/global_delta_min5_cell_qc.R

human-singler:
	conda run -n $(R_ENV) Rscript scripts/signatures/human/run_singler_annotation.R

human-monocytes:
	conda run -n $(R_ENV) Rscript scripts/signatures/human/score_monocyte_classical_nonclassical_signatures.R

human-cluster2-isg:
	conda run -n $(R_ENV) Rscript scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R

pseudobulk:
	conda run -n $(R_ENV) Rscript scripts/pseudobulk/run_cluster_pseudobulk_dge.R
