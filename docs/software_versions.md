# Software Versions And Algorithm References

This table summarizes the main software and algorithms used by the repository. R package versions are taken from `envs/renv.lock`; Cell Ranger is documented by the object-building scripts and environment notes.

| Component | Version used | Repository use | Verified citation/reference note |
| --- | ---: | --- | --- |
| Cell Ranger | 9.0.1 | 10x FASTQ processing with `cellranger multi` for GEX, Feature Barcode/hashtag, and V(D)J libraries | Software/key-resource entry: 10x Genomics Cell Ranger v9.0.1 documentation/release notes; no Bowtie2/SAMtools citation needed for this workflow |
| R | 4.5.1 | Main R runtime | R Core Team |
| Bioconductor | 3.21 | Bioconductor package release | Bioconductor project |
| Seurat | 5.3.0 | Seurat object construction, normalization, dimensional reduction, clustering, plotting, module scores | Hao et al., Nature Biotechnology 2024, doi:10.1038/s41587-023-01767-y |
| SeuratObject | 5.1.0 | Seurat object infrastructure | Satija Lab SeuratObject |
| sctransform | 0.4.2 | SCTransform normalization for the mouse analysis object | Hafemeister and Satija, Genome Biology 2019, doi:10.1186/s13059-019-1874-1; Choudhary and Satija, Genome Biology 2022 |
| SingleR | 2.10.0 | Human cluster annotation against reference data | Aran et al., Nature Immunology 2019, doi:10.1038/s41590-018-0276-y |
| celldex | 1.18.0 | Human Primary Cell Atlas reference used by SingleR | `HumanPrimaryCellAtlasData()` uses HPCA from Mabbott et al., BMC Genomics 2013, processed as described in Aran et al., Nature Immunology 2019 |
| CellChat | 2.2.0 | Cell-cell communication analysis on SCT-normalized expression | CellChat v2 guidance: cite Jin et al., Nature Protocols 2024, doi:10.1038/s41596-024-01045-4; CellChat v1 paper: Jin et al., Nature Communications 2021, doi:10.1038/s41467-021-21246-9 |
| edgeR | 4.6.2 | Sample-level pseudobulk differential expression | Chen et al., Nucleic Acids Research 2025, doi:10.1093/nar/gkaf018 |
| monocle3 | 1.4.26 | Trajectory analysis | Trapnell et al., Nature Biotechnology 2014; Qiu et al., Nature Methods 2017; Cao et al., Nature 2019 |
| orthogene | 1.14.01 | Mouse-human ortholog translation helpers | Schilder and Skene, Bioconductor 2022, doi:10.18129/B9.bioc.orthogene |
| UCell | 2.12.0 | Optional confirmation scoring for selected signatures | Andreatta and Carmona, Computational and Structural Biotechnology Journal 2021, doi:10.1016/j.csbj.2021.06.043 |
| EnhancedVolcano | 1.26.0 | Volcano plot generation | Blighe, Rana, and Lewis, Bioconductor 2025, doi:10.18129/B9.bioc.EnhancedVolcano |
| pheatmap | 1.0.13 | Heatmap generation | Kolde, CRAN 2025, doi:10.32614/CRAN.package.pheatmap |
| ggplot2 | 3.5.2 | Main plotting system | Wickham, ggplot2 |
| ggrepel | 0.9.6 | Volcano plot labels | CRAN `ggrepel` |
| patchwork | 1.3.1 | Multi-panel figure assembly | CRAN `patchwork` |
| dplyr | 1.1.4 | Data manipulation | tidyverse/dplyr |
| tidyr | 1.3.1 | Table reshaping | tidyverse/tidyr |
| readxl | 1.4.5 | Reading edgeR Excel workbooks for selected volcano panels | CRAN `readxl` |
| openxlsx | 4.2.8 | Excel workbook outputs | CRAN `openxlsx` |

Citation checks were performed with `utils::citation()` from the project `renv` library where available. For Cell Ranger, use the 10x Genomics software documentation/release notes as a key-resource/software citation rather than Bowtie2 or SAMtools. For CellChat v2.2.0, package metadata returns a manual citation, but the project README for CellChat v2 asks users of CellChat versions `>= 1.5` to cite the 2024 Nature Protocols paper.

Use `renv::restore(lockfile = "envs/renv.lock")` to recreate the recorded publication package set. A local project library can be newer than the lockfile if packages were installed while `RENV_CONFIG_SYNCHRONIZED_CHECK=false`; the lockfile remains the source of truth for the repository.

## Methods Constants Checked In Scripts

| Analysis stage | Script | Values currently encoded |
| --- | --- | --- |
| Mouse object QC/clustering | `scripts/build_object/build_analysis_seurat_object.R` | `nFeature_RNA > 200`, `nFeature_RNA < 6000`, `percent.mt < 10`, HTODemux positive quantile `0.99`, PCA/UMAP/neighborhood dimensions `1:30`, clustering resolution `1` |
| Human GSE221553 object QC/clustering | `scripts/build_object/build_human_gse221553_seurat_object.R` | `nFeature_RNA > 200`, `nFeature_RNA < 7000`, `percent.mt < 20`, 3000 variable features, PCA/UMAP dimensions `1:30`, clustering resolution `0.5` |
| Human metadata harmonization | `scripts/build_object/build_human_gse221553_seurat_object.R` | GEO sample IDs are converted to `patient_simple`, `timepoint_simple`, and `response_simple`; responders are `P7`, `P8`, `P13`; non-responders are `P1`, `P5`, `P6`, `P10` |
| Human RPCA integration | `scripts/build_object/integrate_human_rpca_by_patient.R` | Seurat RPCA integration by patient, 3000 integration features, dimensions `1:30`, clustering resolution `0.5`; default Figure 7 input is `data/human_seurat_rpca.rds` |
| Biotin thresholds | `scripts/qc/biotin_thresholds_and_status.R` | Untreated cluster-specific `Biotin-TotalSeqC` Q3 threshold; threshold requires at least 5 untreated cells per cluster |
| Biotin delta sensitivity | `scripts/qc/biotin_delta_threshold_sensitivity.R` | Sample/hashtag-level Biotin medians; at least 20 cells per sample-cluster; practical delta guide `0.1` |
| Cluster 7 heatmap | `scripts/visualization/cluster7_marker_heatmaps.R` | Candidate markers from cluster 7 versus clusters 9/8/1/13/7; selected manuscript output is RNA assay, sample-level, filtered genes |
| Human cluster 2 ISG | `scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R` | Final 47-gene high-confidence human ISG ortholog set; primary comparison is responder vs non-responder across all timepoints in cluster 2 |
