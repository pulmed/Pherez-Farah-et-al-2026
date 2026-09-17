#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/human/score_monocyte_classical_nonclassical_signatures.R
# Original file: HUMAN MONO CLASSICAL VS NON CLASSICAL.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Score classical and non-classical monocyte signatures in the human Seurat object.
# Inputs: Human Seurat object with RNA data, UMAP, and cluster metadata.
# Outputs: Feature plots, violin plots, heatmap, gene-presence tables, and optional scored RDS.
# Assay/layer input: RNA data.
# Dependencies: Seurat, SeuratObject, Matrix, dplyr, tidyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Scores are mean normalized expression, matching the uploaded analysis script.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

source("scripts/utils/seurat_io.R")
source("scripts/utils/signature_helpers.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("HUMAN_SEURAT_RDS", unset = "data/human_seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "human", "monocyte_classical_nonclassical"))

assay_use <- Sys.getenv("ASSAY", unset = "RNA")
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "seurat_clusters")
save_scored_rds <- identical(tolower(Sys.getenv("SAVE_SCORED_RDS", unset = "false")), "true")

classical_mono <- c("CD14", "S100A8", "S100A9", "VCAN", "CCR2", "FCN1")
nonclassical_mono <- c("FCGR3A", "CX3CR1", "CDKN1C", "MS4A7")

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(rds_path)) {
  message("SMOKE_TEST=1: human Seurat object not found; skipping human monocyte signature scoring.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# LOAD AND SCORE
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_assays(obj, assay_use)
require_metadata(obj, cluster_col)

DefaultAssay(obj) <- assay_use
obj <- coerce_metadata_character(obj, cluster_col)
obj <- ensure_rna_data_layer(obj, assay = assay_use)

if (length(SeuratObject::Layers(obj[[assay_use]])) > 1) {
  obj[[assay_use]] <- JoinLayers(obj[[assay_use]])
}

expr <- GetAssayData(obj, assay = assay_use, layer = "data")
features <- rownames(expr)

classical_present <- present_signature_genes(classical_mono, features, "ClassicalMono")
nonclassical_present <- present_signature_genes(nonclassical_mono, features, "NonClassicalMono")

gene_check <- bind_rows(
  gene_presence_table(classical_mono, features, "ClassicalMono"),
  gene_presence_table(nonclassical_mono, features, "NonClassicalMono")
)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: classical monocyte genes present: ", length(classical_present), "/", length(classical_mono))
  message("SMOKE_TEST=1: non-classical monocyte genes present: ", length(nonclassical_present), "/", length(nonclassical_mono))
  quit(save = "no", status = 0)
}

write.csv(
  gene_check,
  file.path(outdir, "Monocyte_Classical_vs_NonClassical_gene_presence.csv"),
  row.names = FALSE
)

openxlsx::write.xlsx(
  list(Gene_presence = gene_check),
  file.path(outdir, "Monocyte_Classical_vs_NonClassical_gene_presence.xlsx"),
  overwrite = TRUE
)

obj$ClassicalMono <- Matrix::colMeans(expr[classical_present, , drop = FALSE])
obj$NonClassicalMono <- Matrix::colMeans(expr[nonclassical_present, , drop = FALSE])

# ------------------------------------------------------------------------------
# FEATURE AND VIOLIN PLOTS
# ------------------------------------------------------------------------------
p_feature <- FeaturePlot(
  obj,
  features = c("ClassicalMono", "NonClassicalMono"),
  order = TRUE,
  ncol = 2
)

ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_FeaturePlot.pdf"), p_feature, width = 12, height = 6)
ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_FeaturePlot.png"), p_feature, width = 12, height = 6, dpi = 300)

p_violin <- VlnPlot(
  obj,
  features = c("ClassicalMono", "NonClassicalMono"),
  group.by = cluster_col,
  pt.size = 0,
  ncol = 1
)

ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_ViolinPlot.pdf"), p_violin, width = 12, height = 10)
ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_ViolinPlot.png"), p_violin, width = 12, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# MARKER HEATMAP
# ------------------------------------------------------------------------------
all_mono_genes <- c(classical_present, nonclassical_present)
md <- obj@meta.data
cluster_id <- as.character(md[[cluster_col]])
cluster_levels <- unique(cluster_id)

if (all(grepl("^[0-9]+$", cluster_levels))) {
  cluster_levels <- as.character(sort(as.numeric(cluster_levels)))
}

avg_expr <- lapply(cluster_levels, function(cluster) {
  cells <- rownames(md)[cluster_id == cluster]
  Matrix::rowMeans(expr[all_mono_genes, cells, drop = FALSE])
}) %>%
  do.call(cbind, .)

rownames(avg_expr) <- all_mono_genes
colnames(avg_expr) <- cluster_levels

scaled_expr <- t(scale(t(as.matrix(avg_expr))))
scaled_expr[is.na(scaled_expr)] <- 0

heat_df <- as.data.frame(as.table(scaled_expr))
colnames(heat_df) <- c("gene", "cluster", "zscore")
heat_df <- heat_df %>%
  mutate(
    gene = factor(gene, levels = rev(all_mono_genes)),
    cluster = factor(cluster, levels = cluster_levels)
  )

p_heatmap <- ggplot(heat_df, aes(x = cluster, y = gene, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, name = "Z-score") +
  labs(title = "Classical vs non-classical monocyte markers", x = "Cluster", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_Heatmap.pdf"), p_heatmap, width = 12, height = 7)
ggsave(file.path(outdir, "Monocyte_Classical_vs_NonClassical_Heatmap.png"), p_heatmap, width = 12, height = 7, dpi = 300)

if (save_scored_rds) {
  saveRDS(obj, file.path(outdir, "Human_seurat_with_monocyte_signature_scores.rds"))
}

message("Done. Outputs written to: ", outdir)
