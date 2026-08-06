# Unraveling Crosstalk between T-cells and TME

This repository contains the analysis scripts and computational environment files used in the manuscript **"Unraveling Crosstalk between T-cells and TME."**

## Authors

* **Alfredo Pherez-Farah** - ORCID: 0000-0003-2213-3405
* **Willem de Koning** - ORCID: 0000-0002-4594-8423

## Repository Structure

Scripts are organized by analysis stage instead of by one global numeric sequence. The only linear dependency is object construction; downstream modules consume the same canonical Seurat object and write their own outputs.

```text
Pherez-Farah-et-al-2026/
|-- .Rprofile
|-- .gitignore
|-- README.md
|-- CITATION.cff
|-- LICENSE
|-- Makefile
|-- docs/
|   |-- original_output_audit.md
|   `-- output_manifest.md
|-- envs/
|   |-- README.md
|   |-- environment.yml
|   |-- contam_scan_env.yml
|   `-- renv.lock
|-- renv/
|   |-- .gitignore
|   `-- activate.R
`-- scripts/
    |-- validate_repository_inputs.R
    |-- run_smoke_tests.sh
    |-- utils/
    |   `-- seurat_io.R
    |-- build_object/
    |   |-- cellranger_multi.sh
    |   |-- import_cellranger_multi_to_seurat.R
    |   |-- build_analysis_seurat_object.R
    |   `-- ensure_analysis_layers.R
    |-- qc/
    |   |-- general_metrics_overview.R
    |   |-- biotin_thresholds_and_status.R
    |   `-- scan_contaminants_per_run.sh
    |-- trajectory/
    |   `-- run_trajectory_analysis.R
    |-- cellchat/
    |   |-- run_biotin_comparisons.R
    |   |-- plot_feature_panels_full.R
    |   `-- plot_feature_panels_va2_vb5_ignored.R
    |-- visualization/
    |   |-- cluster7_marker_heatmaps.R
    |   |-- biotin_cellchat_feature_panels.R
    |   `-- biotin_violins_by_cluster.R
    |-- signatures/
    |   |-- score_isg_module_cell_and_sample.R
    |   |-- interrogate_signatures_cell_level.R
    |   |-- compare_m1_m2_signatures_with_violins.R
    |   `-- confirm_signature_trends_sample_level.R
    `-- pseudobulk/
        `-- run_cluster_pseudobulk_dge.R
```

The local `data/`, `output/`, and `original_scripts/` directories are ignored by git. The canonical Seurat RDS is not distributed through GitHub; sequencing/count data should be obtained from GEO and processed locally, or the author-provided object should be placed locally as `data/seurat.rds`. `original_scripts/` is a local provenance archive containing unmodified scripts and result artifacts used during repository cleanup.

## Reviewer Quick Start

The downstream manuscript analyses require a local canonical Seurat object. This RDS is not included in GitHub. Build it from the GEO-deposited counts with the object-building scripts, or place an author-provided annotated object at:

```text
data/seurat.rds
```

The scripts also accept an external object path via `SEURAT_RDS`.

```bash
conda env create -f envs/environment.yml
conda activate pherez-farah-r_env
scripts/run_smoke_tests.sh
Rscript scripts/visualization/biotin_cellchat_feature_panels.R
```

Example with an external object:

```bash
SEURAT_RDS=../Results/R_Output/seurat_final_clean4.rds Rscript scripts/visualization/biotin_cellchat_feature_panels.R
```

## Canonical Seurat Object Contract

The repository is structured around one analysis-ready Seurat object. The object is a local derived input, not a tracked GitHub file. Build it once from GEO counts or place it locally, then use it everywhere:

```text
data/seurat.rds
```

This object corresponds to the original analysis object `../Results/R_Output/seurat_final_clean4.rds`. It is not committed to this repository; public count data are expected to be distributed through GEO. It should be treated as a primary manuscript input, because the final cluster annotation was produced iteratively and is not fully reproduced by the preprocessing scripts.

Required assays and layers:

| Assay | Required layer/slot | Used for |
| --- | --- | --- |
| `RNA` | `counts` | pseudobulk edgeR and count-based summaries |
| `RNA` | `data` | log-normalized expression, feature plots, module scores, heatmaps |
| `SCT` | `data` | SCT marker tests, CellChat expression matrix, some heatmap outputs |
| `ADT` | `counts` and `data` | antibody features, CD8/Biotin gates, biotin thresholding |
| `HTO` | `counts` and/or normalized data when available | hashtag demultiplexing during object construction |

Required reductions and metadata:

| Field | Requirement |
| --- | --- |
| `pca`, `umap` | Required by QC plots, trajectory setup, and feature panels |
| `final_clusters` | Manuscript cluster labels; values must match the author annotation |
| `condition` | Treatment group, expected values `treated` and `untreated` |
| `hash.ID` | Sample identifier for sample-level summaries and pseudobulk designs |

Changing cluster IDs, reclustering, merging, or renaming clusters will change the interpretation of downstream outputs.

## Building the Analysis Object

Raw sequencing data were processed with Cell Ranger 9.0.1. The object-building stage is the only part of the repository that should be read as a sequential pipeline:

| Step | Script | Input | Output |
| --- | --- | --- | --- |
| Cell Ranger multi | `scripts/build_object/cellranger_multi.sh` | FASTQs, 10x references, feature reference CSV | Cell Ranger `*_multi` outputs |
| Import Cell Ranger output | `scripts/build_object/import_cellranger_multi_to_seurat.R` | `CELLRANGER_MULTI_DIR`, default `output/cellranger_multi` | merged multimodal Seurat object |
| Build canonical object | `scripts/build_object/build_analysis_seurat_object.R` | `MERGED_SEURAT_RDS`, default `data/merged_cellranger_multi_seurat.rds` | `SEURAT_RDS`, default `data/seurat.rds` |
| Ensure canonical layers | `scripts/build_object/ensure_analysis_layers.R` | `SEURAT_RDS`, default `data/seurat.rds` | updated object with required normalized layers |

The final manuscript object contains manual/iterative annotation not fully captured by these build scripts. For manuscript reproduction, use the provided/supplied annotated object as `data/seurat.rds`.

## Analysis Modules And Inputs

Each downstream script reads `SEURAT_RDS` by default and writes under `OUTPUT_DIR` or `output/`.

| Analysis | Script | Main input assay/layer | Main outputs |
| --- | --- | --- | --- |
| General metrics/QC | `scripts/qc/general_metrics_overview.R` | `RNA:data`, `SCT:data`, UMAP, metadata | overview plots, cluster composition, marker summaries, selected heatmaps |
| Biotin thresholds/status | `scripts/qc/biotin_thresholds_and_status.R` | `ADT:data` for Biotin thresholding; `RNA:counts` for treated-only edgeR | biotin thresholds, status summaries, legacy Panel 1G delta plots, treated-only pseudobulk DGE |
| Contamination scan | `scripts/qc/scan_contaminants_per_run.sh` | Cell Ranger BAMs and `data/contaminants.fa` | contaminant read/cell-barcode summaries |
| Trajectory | `scripts/trajectory/run_trajectory_analysis.R` | UMAP, `RNA:data`, `final_clusters` | monocle trajectory plots and root-node summaries |
| CellChat comparisons | `scripts/cellchat/run_biotin_comparisons.R` | `SCT:data`, `ADT:data`, cluster/condition metadata | CellChat objects, UMAP grouping plots, communication outputs |
| CellChat feature panels | `scripts/cellchat/plot_feature_panels_full.R` | `RNA:data`, `ADT:data`, UMAP | full comparison feature panels |
| VA2/VB5-ignored feature panels | `scripts/cellchat/plot_feature_panels_va2_vb5_ignored.R` | `RNA:data`, `ADT:data`, UMAP | alternate feature panels |
| Cluster 7 marker heatmaps | `scripts/visualization/cluster7_marker_heatmaps.R` | `RNA:data` and `SCT:data` | cluster 7 marker workbooks and cell/sample heatmaps |
| Biotin CellChat feature panels | `scripts/visualization/biotin_cellchat_feature_panels.R` | `RNA:data`, `ADT:data`, UMAP | biotin/CD8/myeloid feature panels |
| Biotin violins | `scripts/visualization/biotin_violins_by_cluster.R` | `ADT:data`, cluster/condition metadata | per-cluster Biotin violin PDFs |
| ISG module scores | `scripts/signatures/score_isg_module_cell_and_sample.R` | `RNA:data` | cell- and sample-level ISG workbooks/plots |
| Cell-level signatures | `scripts/signatures/interrogate_signatures_cell_level.R` | `RNA:data` | signature ridgeplots, violin plots, selected heatmaps |
| M1/M2 signature comparison | `scripts/signatures/compare_m1_m2_signatures_with_violins.R` | configured expression assay, defaults from object | M1/M2 violin plot and stats workbook |
| Sample-level signature trends | `scripts/signatures/confirm_signature_trends_sample_level.R` | `RNA:data` | sample-level signature heatmaps, ridgeplots, summaries |
| Cluster pseudobulk DGE | `scripts/pseudobulk/run_cluster_pseudobulk_dge.R` | `RNA:counts`, `hash.ID`, `final_clusters` | pairwise edgeR CSV/XLSX tables and volcano PDFs |

## Figure And Result Provenance

| Panel/result | Source |
| --- | --- |
| Panel 1E | `scripts/visualization/cluster7_marker_heatmaps.R`; outputs include `cluster7_top50_RNA.xlsx` and `cluster7_top50_SCT.xlsx` |
| Panel 1F | `scripts/cellchat/run_biotin_comparisons.R`, `scripts/cellchat/plot_feature_panels_full.R`, `scripts/cellchat/plot_feature_panels_va2_vb5_ignored.R`, and `scripts/visualization/biotin_cellchat_feature_panels.R` |
| Panel 1G | `scripts/qc/biotin_thresholds_and_status.R`; outputs include Biotin status summaries and legacy `Delta_Barplot_*` plot names |
| ISG outputs | `scripts/signatures/score_isg_module_cell_and_sample.R` |
| M1/M2 outputs | `scripts/signatures/compare_m1_m2_signatures_with_violins.R` and related signature scripts |
| Pseudobulk DGE tables | `scripts/pseudobulk/run_cluster_pseudobulk_dge.R` |
| Contamination scan/QC | `scripts/qc/scan_contaminants_per_run.sh` |

A detailed check against the unmodified files in `original_scripts/` is recorded in `docs/original_output_audit.md`. A script-to-output map is available in `docs/output_manifest.md`.

## Environment Setup

Environment setup is documented in `envs/README.md`.

In brief:

* R-based analysis uses the main conda environment plus `renv`.
* Contamination scanning uses the separate conda environment defined in `envs/contam_scan_env.yml`.
* Cell Ranger preprocessing requires a working Cell Ranger 9.0.1 installation.

## Results Policy

Results are not stored in this repository. Scripts write to ignored local output directories, usually under `output/`. The supplied `original_scripts/` directory is a local provenance archive and should not be committed.
