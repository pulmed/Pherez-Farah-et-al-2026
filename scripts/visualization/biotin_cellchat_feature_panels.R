#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/visualization/biotin_cellchat_feature_panels.R
# Original file: 20260221 FEATURE PLOT BIOTIN CELL CHAT.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Generate biotin, CD8, and myeloid CellChat feature panels.
# Inputs: Canonical analysis-ready Seurat object with RNA, ADT, UMAP, and cluster metadata.
# Outputs: Biotin CellChat feature panel PDFs.
# Assay/layer input: RNA data and ADT data with UMAP reduction.
# Dependencies: Seurat, ggplot2, dplyr, patchwork, scales.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Uses SEURAT_RDS and OUTPUT_DIR; defaults are data/seurat.rds and output/.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0) CONFIGURATION + LOAD OBJECT
# ------------------------------------------------------------------------------
seurat_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
seurat_rdata <- Sys.getenv("SEURAT_RDATA", unset = file.path(dirname(seurat_rds), "seurat_final_clean4.RData"))
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")

if (file.exists(seurat_rds)) {
  obj <- readRDS(seurat_rds)
} else if (file.exists(seurat_rdata)) {
  load(seurat_rdata)
  if (!exists("obj")) {
    if (exists("seurat_final_clean4")) {
      obj <- seurat_final_clean4
    } else {
      stop("Loaded seurat_final_clean4.RData but neither 'obj' nor 'seurat_final_clean4' exists in the environment.")
    }
  }
} else {
  stop("Could not find Seurat input. Set SEURAT_RDS or place seurat_final_clean4.rds at ", seurat_rds)
}

stopifnot(inherits(obj, "Seurat"))
stopifnot("final_clusters" %in% colnames(obj@meta.data))
stopifnot("condition" %in% colnames(obj@meta.data))
stopifnot("umap" %in% names(obj@reductions))

obj$final_clusters <- as.character(obj$final_clusters)

# ------------------------------------------------------------------------------
# Output directory (as requested)
# ------------------------------------------------------------------------------
outdir <- file.path(output_root, "20260221 FEATURE PLOT WITH BIOTIN")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1) COMPUTE biotin_status IF MISSING
# ------------------------------------------------------------------------------
BIOTIN_ADT_FEATURE <- "Biotin-TotalSeqC"
MIN_UNTREATED_PER_CLUSTER <- 5
Q_BIOTIN <- 0.75

if (!("biotin_status" %in% colnames(obj@meta.data)) || all(is.na(obj$biotin_status))) {

  message("biotin_status missing/empty -> computing from cluster-specific Q3 of untreated...")

  DefaultAssay(obj) <- "ADT"

  obj <- NormalizeData(
    obj,
    normalization.method = "CLR",
    margin = 2,
    verbose = FALSE
  )

  adt <- GetAssayData(obj, assay = "ADT", slot = "data")
  if (!(BIOTIN_ADT_FEATURE %in% rownames(adt))) {
    stop("Missing ADT feature: ", BIOTIN_ADT_FEATURE)
  }

  meta <- obj@meta.data
  meta$final_clusters <- as.character(meta$final_clusters)

  clusters <- sort(unique(meta$final_clusters))
  thresholds <- setNames(rep(NA_real_, length(clusters)), clusters)

  # thresholds from untreated cells only
  for (cl in clusters) {
    cells_u <- rownames(meta)[meta$final_clusters == cl & meta$condition == "untreated"]
    if (length(cells_u) >= MIN_UNTREATED_PER_CLUSTER) {
      thresholds[[cl]] <- as.numeric(quantile(adt[BIOTIN_ADT_FEATURE, cells_u], Q_BIOTIN, na.rm = TRUE))
    }
  }

  # assign status for all cells in each cluster (where threshold exists)
  biotin_status <- rep(NA_character_, ncol(obj))
  names(biotin_status) <- colnames(obj)

  for (cl in names(thresholds)) {
    thr <- thresholds[[cl]]
    if (is.na(thr)) next
    cells_cl <- rownames(meta)[meta$final_clusters == cl]
    vals <- adt[BIOTIN_ADT_FEATURE, cells_cl]
    biotin_status[cells_cl] <- ifelse(vals > thr, "positive", "negative")
  }

  obj$biotin_status <- factor(biotin_status, levels = c("negative", "positive"))
  DefaultAssay(obj) <- "RNA"

  message("biotin_status computed. Counts:")
  print(table(obj$biotin_status, useNA = "ifany"))
} else {
  message("biotin_status found; using existing column.")
  obj$biotin_status <- factor(as.character(obj$biotin_status), levels = c("negative", "positive"))
}

# ------------------------------------------------------------------------------
# 2) SETTINGS: CLUSTERS + COLORS + HELPERS
# ------------------------------------------------------------------------------
tcell_clusters <- c("6", "10", "11")
myeloid_clusters <- as.character(c(0,1,2,3,4,5,7,8,9,13,14,16,17))

COL_OTHER <- "grey75"
COL_BLUE <- "#859EE1"   # requested blue
COL_RED <- "#ea6d6d"
COL_YELLOW <- "#FFD84D"   # overlay for treated myeloid biotin+

get_cluster_centers <- function(objx, reduction = "umap", cluster_col = "final_clusters") {
  emb <- Embeddings(objx, reduction = reduction)
  data.frame(
    x = emb[, 1],
    y = emb[, 2],
    cluster = as.character(objx[[cluster_col]][, 1]),
    row.names = rownames(emb)
  ) %>%
    group_by(cluster) %>%
    summarise(x = median(x), y = median(y), .groups = "drop")
}
centers <- get_cluster_centers(obj)

meta <- obj@meta.data
is_treated <- meta$condition == "treated"
Biotin_pos <- meta$biotin_status == "positive"

treated_myeloid_biotinpos_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% myeloid_clusters &
    Biotin_pos
]

# ------------------------------------------------------------------------------
# 3) PLOT A: Treated Myeloid biotin+ in BLUE + cluster labels
# ------------------------------------------------------------------------------
obj$plot_group_biotin_myel <- "Other"
obj$plot_group_biotin_myel[colnames(obj) %in% treated_myeloid_biotinpos_cells] <- "Treated Myeloid biotin+"
obj$plot_group_biotin_myel <- factor(obj$plot_group_biotin_myel, levels = c("Other", "Treated Myeloid biotin+"))

p1 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "plot_group_biotin_myel",
  pt.size = 0.25,
  cols = c("Other" = COL_OTHER, "Treated Myeloid biotin+" = COL_BLUE)
) +
  ggtitle("Treated myeloid clusters: Biotin+ (threshold = cluster Q3 of untreated)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  geom_text(
    data = centers,
    aes(x = x, y = y, label = cluster),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  )

pdf(file.path(outdir, "UMAP_TreatedMyeloid_BiotinPos_Blue.pdf"), width = 7.5, height = 6.5)
print(p1)
dev.off()

# ------------------------------------------------------------------------------
# 4) PLOT B: Comp1B-style + YELLOW overlay for treated myeloid biotin+
# ------------------------------------------------------------------------------

# OT-I definition (as in your final pipeline):
# OT-I = treated & T clusters & CD8_gate & biotin+ & (Va2+Vb5+)

VA2_ADT_FEATURE <- "Va2-TotalSeqC"
VB5_ADT_FEATURE <- "Vb5-TotalSeqC"
CD8A_ADT_FEATURE <- "CD8a-TotalSeqC"

# RNA CD8 (SCT preferred; robust)
get_gene_vec <- function(objx, gene, prefer_assays = c("SCT", "RNA")) {
  assay_use <- NULL
  for (a in prefer_assays) {
    if (a %in% names(objx@assays) && gene %in% rownames(objx[[a]])) { assay_use <- a; break }
  }
  if (is.null(assay_use)) return(setNames(rep(0, ncol(objx)), colnames(objx)))

  m <- tryCatch(GetAssayData(objx, assay = assay_use, slot = "data"), error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0 || !(gene %in% rownames(m))) {
    return(setNames(rep(0, ncol(objx)), colnames(objx)))
  }
  v <- as.numeric(m[gene, , drop = TRUE])
  names(v) <- colnames(objx)
  v
}

Cd8a_expr <- get_gene_vec(obj, "Cd8a")
Cd8b1_expr <- get_gene_vec(obj, "Cd8b1")
CD8_rna_pos <- (Cd8a_expr > 0) | (Cd8b1_expr > 0)

DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
adt <- GetAssayData(obj, assay = "ADT", slot = "data")

need <- c(VA2_ADT_FEATURE, VB5_ADT_FEATURE, CD8A_ADT_FEATURE)
missing <- setdiff(need, rownames(adt))
if (length(missing) > 0) stop("Missing ADT features: ", paste(missing, collapse = ", "))

Va2_pos <- adt[VA2_ADT_FEATURE, ]  > 0
Vb5_pos <- adt[VB5_ADT_FEATURE, ]  > 0
CD8a_pos <- adt[CD8A_ADT_FEATURE, ] > 0
CD8_gate <- CD8a_pos | CD8_rna_pos

meta <- obj@meta.data
meta$final_clusters <- as.character(meta$final_clusters)
is_treated <- meta$condition == "treated"
Biotin_pos <- meta$biotin_status == "positive"

OTI_comp1B_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% tcell_clusters &
    CD8_gate &
    Biotin_pos &
    (Va2_pos & Vb5_pos)
]

Myeloid_treated_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% myeloid_clusters
]

# Priority group assignment:
#  - OT-I (red)
#  - Myeloid biotin+ treated (yellow)
#  - Myeloid treated (blue)
#  - Other (grey)
grp <- rep("Other", ncol(obj))
names(grp) <- colnames(obj)

grp[colnames(obj) %in% Myeloid_treated_cells] <- "Myeloid"
grp[colnames(obj) %in% treated_myeloid_biotinpos_cells] <- "Myeloid biotin+ (treated)"
grp[colnames(obj) %in% OTI_comp1B_cells] <- "OT-I"

obj$plot_group_comp1B_overlay <- factor(
  grp,
  levels = c("Other", "Myeloid", "Myeloid biotin+ (treated)", "OT-I")
)

DefaultAssay(obj) <- "RNA"

p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "plot_group_comp1B_overlay",
  pt.size = 0.25,
  cols = c(
    "Other" = COL_OTHER,
    "Myeloid" = COL_BLUE,
    "Myeloid biotin+ (treated)" = COL_YELLOW,
    "OT-I" = COL_RED
  )
) +
  ggtitle("Comp1B overlay: treated myeloid biotin+ (yellow) + OT-I (red)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  geom_text(
    data = centers,
    aes(x = x, y = y, label = cluster),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  )

pdf(file.path(outdir, "UMAP_Comp1B_withMyeloidBiotinPos_YellowOverlay.pdf"), width = 7.5, height = 6.5)
print(p2)
dev.off()

message("Done. Saved PDFs in: ", outdir)

