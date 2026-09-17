# Output Manifest

This manifest lists the main analysis scripts, their expected Seurat assay/layer inputs, and their default output locations. Most scripts accept `SEURAT_RDS` and `OUTPUT_DIR` overrides.

| Analysis | Script | Main input | Default output | Manuscript relevance |
| --- | --- | --- | --- | --- |
| Build canonical Seurat object | `scripts/build_object/build_analysis_seurat_object.R` | merged Seurat object; RNA/ADT/HTO counts | local `data/seurat.rds`, `output/build_analysis_seurat_object/` | canonical downstream input |
| Ensure canonical layers | `scripts/build_object/ensure_analysis_layers.R` | local `data/seurat.rds` | updated local `data/seurat.rds` | guarantees required normalized layers |
| Build human GSE221553 object | `scripts/build_object/build_human_gse221553_seurat_object.R` | local or downloaded GSE221553 `*-counts.tsv.gz` tables | local `data/human_seurat.rds`, `output/human/build_gse221553/` | human validation object provenance |
| General metrics | `scripts/qc/general_metrics_overview.R` | `RNA:data`, `SCT:data`, UMAP, metadata | `output/03_general_metrics_overview/` | QC and exploratory summaries |
| Biotin thresholds/status | `scripts/qc/biotin_thresholds_and_status.R` | `ADT:data`, `RNA:counts`, metadata | `output/05_biotin_analysis/` | Panel 1G and treated-only biotin DGE |
| Contamination scan | `scripts/qc/scan_contaminants_per_run.sh` | Cell Ranger BAMs, `data/contaminants.fa` | `output/contam_per_run/` | QC/provenance |
| Trajectory analysis | `scripts/trajectory/run_trajectory_analysis.R` | UMAP, `RNA:data`, `final_clusters` | `output/04_trajectory_analysis/` | trajectory/pseudotime outputs |
| CellChat comparisons | `scripts/cellchat/run_biotin_comparisons.R` | `SCT:data`, `ADT:data`, metadata | `output/06_cellchat_biotin_normalized/` | CellChat comparison outputs |
| Full CellChat feature panels | `scripts/cellchat/plot_feature_panels_full.R` | `RNA:data`, `ADT:data`, UMAP | `output/07_CELLCHAT_FEATURES_*` | CellChat feature visualizations |
| VA2/VB5-ignored panels | `scripts/cellchat/plot_feature_panels_va2_vb5_ignored.R` | `RNA:data`, `ADT:data`, UMAP | `output/08_CELLCHAT_FEATURES_*` | alternate CellChat visualizations |
| Biotin CellChat feature panels | `scripts/visualization/biotin_cellchat_feature_panels.R` | `RNA:data`, `ADT:data`, UMAP | `output/20260221 FEATURE PLOT WITH BIOTIN/` | Panel 1F provenance |
| Cluster 7 marker heatmaps | `scripts/visualization/cluster7_marker_heatmaps.R` | `RNA:data`, `SCT:data`, `hash.ID` | `output/11_cluster7_top50_global_markers_compare_9_8_1_13_7/` | Panel 1E provenance |
| Biotin violins | `scripts/visualization/biotin_violins_by_cluster.R` | `ADT:data`, metadata | `output/20260221 BIOTIN VIOLINS_*` | Biotin visualization provenance |
| ISG module scores | `scripts/signatures/score_isg_module_cell_and_sample.R` | `RNA:data`, metadata | `output/12_20260218_ISG_MODULE/` | ISG signature outputs |
| Cross-presentation signature | `scripts/signatures/score_cross_presentation_signature.R` | `RNA:data`, cluster/condition/sample metadata | `output/17_cross_presentation_signature/` | final cross-presentation signature outputs |
| Mouse-to-human ortholog translation | `scripts/signatures/cross_species/translate_mouse_signatures_to_human.R` | final mouse signature vectors | `output/cross_species/ortholog_translation/` | optional cross-species transparency table |
| Cell-level signatures | `scripts/signatures/interrogate_signatures_cell_level.R` | `RNA:data`, UMAP, metadata | `output/13_*` | signature visualization outputs |
| M1/M2 signature comparison | `scripts/signatures/compare_m1_m2_signatures_with_violins.R` | `RNA:data`, cluster/sample metadata | `output/14_M1_VS_M2/M1_M2_final/` | final M1/M2 signature outputs |
| Sample-level signature trends | `scripts/signatures/confirm_signature_trends_sample_level.R` | `RNA:data`, `hash.ID`, metadata | `output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/` | sample-level signature confirmation |
| Human global delta QC | `scripts/qc/human/global_delta_min5_cell_qc.R` | human object metadata | `output/human/global_delta_min5_cell_qc/` | human T0/T30 cluster-fraction min-cell QC |
| Human monocyte signatures | `scripts/signatures/human/score_monocyte_classical_nonclassical_signatures.R` | human `RNA:data`, UMAP, cluster metadata | `output/human/monocyte_classical_nonclassical/` | human monocyte signature interrogation |
| Mouse-to-human cluster 2 ISG | `scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R` | human `RNA:data`, patient/timepoint/response/cluster metadata | `output/human/cluster2_isg_responder_status/` | human cluster 2 responder/non-responder IFN signature comparison |
| Cluster pseudobulk DGE | `scripts/pseudobulk/run_cluster_pseudobulk_dge.R` | `RNA:counts`, `hash.ID`, `final_clusters` | `output/16_DGE_clusters_RNA_pseudobulk/` | pseudobulk DGE tables |
