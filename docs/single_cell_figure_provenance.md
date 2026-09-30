# Single-Cell Figure Provenance

This document maps the single-cell and CITE-seq figure panels to the repository scripts that generate the underlying analyses. Generated results are written to ignored local `output/` directories and are not committed to GitHub.

## Primary Inputs

| Input | Local default | Notes |
| --- | --- | --- |
| Mouse analysis-ready Seurat object | `data/seurat.rds` | Derived from Cell Ranger outputs and final manual/iterative annotation; not committed to GitHub |
| Human validation Seurat object | `data/human_seurat.rds` | Rebuildable from GEO `GSE221553` count tables with repository scripts; not committed to GitHub |
| Human RPCA object | `data/human_seurat_rpca.rds` | Default `HUMAN_ANALYSIS_RDS` for Figure 7 scripts; built from the human object by RPCA integration by patient |

## Object Construction And Preprocessing

| Step | Script | Main output |
| --- | --- | --- |
| Cell Ranger multi preprocessing | `scripts/build_object/cellranger_multi.sh` | Cell Ranger multi outputs for GEX, ADT/HTO, and V(D)J inputs |
| Mouse Cell Ranger import | `scripts/build_object/import_cellranger_multi_to_seurat.R` | `data/merged_cellranger_multi_seurat.rds` |
| Mouse canonical object construction | `scripts/build_object/build_analysis_seurat_object.R` | `data/seurat.rds` |
| Mouse layer check/repair | `scripts/build_object/ensure_analysis_layers.R` | Updated `data/seurat.rds` if normalized layers are missing |
| Human GEO object construction | `scripts/build_object/build_human_gse221553_seurat_object.R` | `data/human_seurat.rds` |
| Human RPCA integration | `scripts/build_object/integrate_human_rpca_by_patient.R` | `data/human_seurat_rpca.rds` |

## Figure 2 And Related Supplements

| Panel/result | Repository script | Main output/provenance |
| --- | --- | --- |
| Figure 2B-C, mouse UMAPs by cluster and condition | `scripts/qc/general_metrics_overview.R` | `output/03_general_metrics_overview/umap_overview/` |
| Figure 2D, therapy-associated cluster representation | `scripts/qc/general_metrics_overview.R` | `output/03_general_metrics_overview/cluster_condition_composition/` |
| Figure 2E/H and Figure S2, RNA/ADT marker visualization | `scripts/qc/general_metrics_overview.R`; `scripts/cellchat/plot_feature_panels_full.R`; `scripts/visualization/biotin_cellchat_feature_panels.R` | Feature plots, violin plots, and selected marker panels under `output/03_*`, `output/07_*`, and `output/20260221 FEATURE PLOT WITH BIOTIN/` |
| Figure 2F, macrophage pseudobulk DGE volcano plots | `scripts/pseudobulk/run_cluster_pseudobulk_dge.R`; `scripts/visualization/plot_pseudobulk_volcano_panels.R` | edgeR workbooks and volcano PDFs under `output/16_DGE_clusters_RNA_pseudobulk/` and selected panels under `output/15_pseudobulk_volcano_panels/` |
| Figure 2G, LAM signature ridge/violin plots | `scripts/signatures/interrogate_signatures_cell_level.R`; `scripts/signatures/confirm_signature_trends_sample_level.R` | LAM signature outputs under `output/13_*` and sample-level confirmation outputs |
| Figure 2I, mregDC/ISG-DC signature scoring | `scripts/signatures/interrogate_signatures_cell_level.R`; `scripts/signatures/confirm_signature_trends_sample_level.R` | DC signature score plots and summaries under `output/13_*` and sample-level confirmation outputs |
| Figure 2J, cluster 7 heatmap | `scripts/visualization/cluster7_marker_heatmaps.R` | Selected manuscript output: `output/11_cluster7_top50_global_markers_compare_9_8_1_13_7/RNA/HEATMAP_sample_level_FILTERED_RNA.pdf`; workbook: `cluster7_top50_RNA.xlsx` |
| Figure 2K and Figure S3A-B, trajectory/IFN program | `scripts/trajectory/run_trajectory_analysis.R`; `scripts/signatures/score_isg_module_cell_and_sample.R` | Monocle trajectory outputs under `output/04_trajectory_analysis/`; ISG module outputs under `output/12_20260218_ISG_MODULE/` |

## Figure 3

| Panel/result | Repository script | Main output/provenance |
| --- | --- | --- |
| Figure 3A, cluster-specific Biotin thresholds | `scripts/qc/biotin_thresholds_and_status.R` | `output/05_biotin_analysis/thresholds/biotin_thresholds_by_cluster.csv` |
| Figure 3B, sample-level Biotin delta waterfall/bar plot | `scripts/qc/biotin_thresholds_and_status.R` | `output/05_biotin_analysis/plots/Delta_Barplot_SampleLevel_SelectedClusters.*`; source table `output/05_biotin_analysis/summaries/biotin_delta_sample_level.csv` |
| Figure 3B support for 0.1 delta guide | `scripts/qc/biotin_delta_threshold_sensitivity.R` | Noise-floor, bootstrap-CI, and effect-size tables under `output/05_biotin_analysis/delta_threshold_sensitivity/` |
| Figure 3C, OT-I/myeloid population definition feature panels | `scripts/cellchat/plot_feature_panels_full.R`; `scripts/cellchat/plot_feature_panels_va2_vb5_ignored.R`; `scripts/visualization/biotin_cellchat_feature_panels.R` | Feature-panel PDFs under `output/07_*`, `output/08_*`, and `output/20260221 FEATURE PLOT WITH BIOTIN/` |
| Figure 3D, CellChat chord plots | `scripts/cellchat/run_biotin_comparisons.R` | CellChat objects and communication outputs under `output/06_cellchat_biotin_normalized/` |
| Figure 3E, cross-presentation signature | `scripts/signatures/score_cross_presentation_signature.R`; optional confirmation in `scripts/signatures/confirm_ucell_signature_scores.R` | Primary AddModuleScore outputs under `output/17_cross_presentation_signature/`; optional UCell confirmation under `output/18_ucell_signature_confirmation/` |

## Figure 7 And Related Supplements

| Panel/result | Repository script | Main output/provenance |
| --- | --- | --- |
| Figure 7A, human myeloid UMAP and annotation | `scripts/build_object/build_human_gse221553_seurat_object.R`; `scripts/build_object/integrate_human_rpca_by_patient.R`; `scripts/signatures/human/run_singler_annotation.R` | Human object, RPCA object, and SingleR/celldex annotations under `output/human/singler_annotation/` |
| Figure 7B, min-5-cell patient-level abundance delta QC | `scripts/qc/human/global_delta_min5_cell_qc.R` | `output/human/global_delta_min5_cell_qc/GLOBAL_DELTA_MIN5_CELL_COUNT_QC.xlsx` and waterfall tables/plots |
| Figure 7C, mouse-to-human IFN signature in human cluster 2 | `scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R` | `output/human/cluster2_isg_responder_status/`; final 47-gene high-confidence ortholog set |
| Figure 7D-E, human response signatures interrogated in mouse | `scripts/signatures/cross_species/translate_mouse_signatures_to_human.R`; `scripts/signatures/interrogate_signatures_cell_level.R`; `scripts/signatures/confirm_signature_trends_sample_level.R` | Ortholog tables under `output/cross_species/ortholog_translation/` and mouse signature plots/summaries under `output/13_*` and sample-level confirmation outputs |
| Figure S6, human monocyte classical/non-classical signatures | `scripts/signatures/human/score_monocyte_classical_nonclassical_signatures.R` | `output/human/monocyte_classical_nonclassical/` |

## Notes On Non-Single-Cell Figures

Figures 1, 4, 5, and 6 are primarily flow-cytometry and experimental perturbation figures. Their statistical tests and replicate definitions are described in the manuscript legends; they are not reconstructed from the single-cell Seurat object by the repository scripts.

## Manuscript-Relevant Analysis Choices

| Choice | Repository location |
| --- | --- |
| Mouse QC/clustering constants | `docs/software_versions.md`; `scripts/build_object/build_analysis_seurat_object.R` |
| Human QC/RPCA constants | `docs/software_versions.md`; `scripts/build_object/build_human_gse221553_seurat_object.R`; `scripts/build_object/integrate_human_rpca_by_patient.R` |
| Raw/log-normalized/SCT/ADT layer use by analysis | README section "Canonical Seurat Object Contract" and "Analysis Modules And Inputs" |
| Final 47-gene human IFN/ISG set | `scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R` |
| Optional UCell confirmation, not primary scoring | `scripts/signatures/confirm_ucell_signature_scores.R` |
