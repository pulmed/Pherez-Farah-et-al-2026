#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/human/run_singler_annotation.R
# Original file: 20260929 SINGLER HUMAN.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Annotate human clusters using SingleR and the celldex Human Primary Cell Atlas reference.
# Inputs: Human Seurat object with RNA data and cluster metadata.
# Outputs: SingleR prediction objects, cluster annotation tables, plots, and annotated RDS.
# Assay/layer input: RNA data.
# Dependencies: Seurat, SingleR, celldex, ggplot2, dplyr, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - SingleR is run at cluster level to reduce memory use.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
})

source("scripts/utils/seurat_io.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv(
  "HUMAN_ANALYSIS_RDS",
  unset = Sys.getenv("HUMAN_RPCA_RDS", unset = "data/human_seurat_rpca.rds")
)
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "human", "singler_annotation"))
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "seurat_clusters")
save_annotated_rds <- identical(tolower(Sys.getenv("SAVE_ANNOTATED_RDS", unset = "false")), "true")

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(rds_path)) {
  message("SMOKE_TEST=1: human Seurat object not found; skipping SingleR annotation.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# LOAD AND VALIDATE
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_assays(obj, "RNA")
require_metadata(obj, cluster_col)

DefaultAssay(obj) <- "RNA"
if (length(SeuratObject::Layers(obj[["RNA"]])) > 1) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
obj <- ensure_rna_data_layer(obj, assay = "RNA")

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: human SingleR cluster column: ", cluster_col)
  message("SMOKE_TEST=1: cells: ", ncol(obj))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# SINGLER ANNOTATION
# ------------------------------------------------------------------------------
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
ref <- celldex::HumanPrimaryCellAtlasData()

clusters <- as.character(obj@meta.data[[cluster_col]])
names(clusters) <- colnames(obj)

pred_broad <- SingleR(
  test = data_mat,
  ref = ref,
  labels = ref$label.main,
  clusters = clusters
)

pred_fine <- SingleR(
  test = data_mat,
  ref = ref,
  labels = ref$label.fine,
  clusters = clusters
)

broad_labels <- pred_broad$pruned.labels
broad_labels[is.na(broad_labels)] <- "Unassigned"
names(broad_labels) <- rownames(pred_broad)

fine_labels <- pred_fine$pruned.labels
fine_labels[is.na(fine_labels)] <- "Unassigned"
names(fine_labels) <- rownames(pred_fine)

obj$SingleR_broad <- unname(broad_labels[clusters])
obj$SingleR_fine <- unname(fine_labels[clusters])

annotation_table <- data.frame(
  cluster = rownames(pred_broad),
  broad_label = pred_broad$labels,
  broad_pruned = pred_broad$pruned.labels,
  fine_label = pred_fine$labels,
  fine_pruned = pred_fine$pruned.labels,
  stringsAsFactors = FALSE
)
annotation_table$n_cells <- as.integer(table(factor(clusters, levels = annotation_table$cluster)))

saveRDS(pred_broad, file.path(outdir, "SingleR_broad_predictions.rds"))
saveRDS(pred_fine, file.path(outdir, "SingleR_fine_predictions.rds"))
write.csv(annotation_table, file.path(outdir, "SingleR_cluster_annotations.csv"), row.names = FALSE)
openxlsx::write.xlsx(
  list(SingleR_cluster_annotations = annotation_table),
  file.path(outdir, "SingleR_cluster_annotations.xlsx"),
  overwrite = TRUE
)

if ("umap" %in% Reductions(obj)) {
  p_broad <- DimPlot(obj, reduction = "umap", group.by = "SingleR_broad", label = TRUE, repel = TRUE) +
    ggtitle("SingleR broad labels")
  p_fine <- DimPlot(obj, reduction = "umap", group.by = "SingleR_fine", label = TRUE, repel = TRUE) +
    ggtitle("SingleR fine labels")
  ggsave(file.path(outdir, "SingleR_broad_umap.pdf"), p_broad, width = 8, height = 6)
  ggsave(file.path(outdir, "SingleR_fine_umap.pdf"), p_fine, width = 8, height = 6)
}

if (save_annotated_rds) {
  saveRDS(obj, file.path(outdir, "human_seurat_with_SingleR_annotations.rds"))
}

message("Done. SingleR outputs written to: ", outdir)
