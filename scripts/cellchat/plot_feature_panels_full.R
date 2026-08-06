#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/cellchat/plot_feature_panels_full.R
# Original file: 20260221 FEATURE PLOTS CELL CHAT.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Generate full CellChat-oriented feature panels for the primary comparison design.
# Inputs: Canonical analysis-ready Seurat object with RNA, ADT, UMAP, and cluster metadata.
# Outputs: Feature panel PDFs for the full CellChat comparison design.
# Assay/layer input: RNA data and ADT data with UMAP reduction.
# Dependencies: Seurat, ggplot2, dplyr, patchwork, scales.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Output paths are rooted under OUTPUT_DIR, default output/.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!exists("obj")) {
  seurat_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
  obj <- readRDS(seurat_rds)
}
stopifnot(exists("obj"))
stopifnot(inherits(obj, "Seurat"))
stopifnot("final_clusters" %in% colnames(obj@meta.data))
stopifnot("condition"      %in% colnames(obj@meta.data))
stopifnot("umap"           %in% names(obj@reductions))

# Ensure cluster column is character (robust)
obj$final_clusters <- as.character(obj$final_clusters)

# ------------------------------------------------------------------------------
# Output folder (new subfolder)
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "07_CELLCHAT_FEATURES_v2_biotinOnlyOTI_CD8RNA_fix")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1) Biotin status (cluster-specific Q3 from UNTREATED) - compute only if missing/empty
# ------------------------------------------------------------------------------
compute_biotin_status_if_missing <- function(seurat_object,
                                             cluster_col = "final_clusters",
                                             condition_col = "condition",
                                             biotin_adt = "Biotin-TotalSeqC",
                                             untreated_label = "untreated",
                                             min_cells = 5,
                                             q = 0.75) {
  if ("biotin_status" %in% colnames(seurat_object@meta.data)) {
    v <- seurat_object$biotin_status
    if (sum(!is.na(v)) > 0) {
      seurat_object$biotin_status <- factor(as.character(v), levels = c("negative", "positive"))
      return(seurat_object)
    }
  }

  DefaultAssay(seurat_object) <- "ADT"
  seurat_object <- NormalizeData(seurat_object, normalization.method = "CLR", margin = 2, verbose = FALSE)

  adt_data <- GetAssayData(seurat_object, assay = "ADT", slot = "data")
  stopifnot(biotin_adt %in% rownames(adt_data))

  get_biotin <- function(objx, cells) {
    GetAssayData(objx, assay = "ADT", slot = "data")[biotin_adt, cells, drop = TRUE]
  }

  meta <- seurat_object@meta.data
  clusters <- sort(unique(as.character(meta[[cluster_col]])))
  thresholds <- setNames(rep(NA_real_, length(clusters)), clusters)

  for (cl in clusters) {
    cells_u <- rownames(meta)[meta[[cluster_col]] == cl & meta[[condition_col]] == untreated_label]
    if (length(cells_u) >= min_cells) {
      vals <- get_biotin(seurat_object, cells_u)
      thresholds[[cl]] <- as.numeric(quantile(vals, q, na.rm = TRUE))
    }
  }

  biotin_status <- rep(NA_character_, ncol(seurat_object))
  names(biotin_status) <- colnames(seurat_object)

  for (cl in names(thresholds)) {
    thr <- thresholds[[cl]]
    if (is.na(thr)) next
    cells_cl <- rownames(meta)[meta[[cluster_col]] == cl]
    if (length(cells_cl) == 0) next
    vals <- get_biotin(seurat_object, cells_cl)
    biotin_status[names(vals)] <- ifelse(vals > thr, "positive", "negative")
  }

  seurat_object$biotin_status <- factor(biotin_status, levels = c("negative", "positive"))
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object
}

obj <- compute_biotin_status_if_missing(obj)

# ------------------------------------------------------------------------------
# 2) Robust gene extractor for RNA/SCT (fixes "Layer data is empty")
# ------------------------------------------------------------------------------
get_gene_vec <- function(objx, gene, prefer_assays = c("SCT", "RNA")) {
  assay_use <- NULL
  for (a in prefer_assays) {
    if (a %in% names(objx@assays) && gene %in% rownames(objx[[a]])) {
      assay_use <- a
      break
    }
  }
  if (is.null(assay_use)) {
    return(setNames(rep(0, ncol(objx)), colnames(objx)))
  }

  pull_mat <- function(slot_or_layer) {
    m <- tryCatch(
      GetAssayData(objx, assay = assay_use, slot = slot_or_layer),
      error = function(e) NULL
    )
    if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
      m <- tryCatch(
        GetAssayData(objx, assay = assay_use, layer = slot_or_layer),
        error = function(e) NULL
      )
    }
    m
  }

  m_data <- pull_mat("data")
  m_counts <- pull_mat("counts")

  m <- m_data
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) m <- m_counts

  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0 || !(gene %in% rownames(m))) {
    return(setNames(rep(0, ncol(objx)), colnames(objx)))
  }

  v <- m[gene, , drop = TRUE]
  v <- as.numeric(v)
  names(v) <- colnames(objx)
  v
}

Cd8a_expr <- get_gene_vec(obj, "Cd8a")
Cd8b1_expr <- get_gene_vec(obj, "Cd8b1")
Cd8b_expr <- get_gene_vec(obj, "Cd8b")

CD8_rna_pos <- (Cd8a_expr > 0) | (Cd8b1_expr > 0) | (Cd8b_expr > 0)

# ------------------------------------------------------------------------------
# 3) ADT gates (CLR) + updated CD8 gate
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
adt <- GetAssayData(obj, assay = "ADT", slot = "data")

needed_adt <- c("Va2-TotalSeqC", "Vb5-TotalSeqC", "CD8a-TotalSeqC")
missing_adt <- setdiff(needed_adt, rownames(adt))
if (length(missing_adt) > 0) stop("Missing ADT features: ", paste(missing_adt, collapse = ", "))

Va2_pos <- adt["Va2-TotalSeqC", ]  > 0
Vb5_pos <- adt["Vb5-TotalSeqC", ]  > 0
CD8a_pos <- adt["CD8a-TotalSeqC", ] > 0

Va2_neg <- !Va2_pos
Vb5_neg <- !Vb5_pos
Va2_xor_Vb5 <- (Va2_pos & !Vb5_pos) | (!Va2_pos & Vb5_pos)

# UPDATED: CD8 gate can be met via ADT or RNA/SCT expression
CD8_gate <- CD8a_pos | CD8_rna_pos

Biotin_pos <- obj$biotin_status == "positive"
Biotin_neg <- obj$biotin_status == "negative"

meta <- obj@meta.data
is_treated <- meta$condition == "treated"
is_untreated <- meta$condition == "untreated"

tcell_clusters <- c("6", "10", "11")  # exclude 15
myeloid_clusters <- as.character(c(0,1,2,3,4,5,7,8,9,13,14,16,17))

# ------------------------------------------------------------------------------
# 4) Six definitions (as used in plots)
# ------------------------------------------------------------------------------
# Comp1B T/OT: OT-I treated, Va2+Vb5+, CD8_gate, AND biotin_pos (ONLY here)
OTI_comp1B_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% tcell_clusters &
    Va2_pos & Vb5_pos &
    CD8_gate &
    Biotin_pos
]

# Comp1B myeloid: treated myeloid/DC, NO biotin filter
Myeloid_treated_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% myeloid_clusters
]

# Comp2 T: Endog CD8 treated, CD8_gate, biotin_neg, and (Va2/Vb5 double-neg OR XOR)
EndogCD8_treated_cells <- colnames(obj)[
  is_treated &
    meta$final_clusters %in% tcell_clusters &
    CD8_gate &
    Biotin_neg &
    ((Va2_neg & Vb5_neg) | Va2_xor_Vb5)
]

# Comp3 T: Endog CD8 untreated, CD8_gate, biotin_neg, and (double-neg OR XOR)
EndogCD8_untreated_cells <- colnames(obj)[
  is_untreated &
    meta$final_clusters %in% tcell_clusters &
    CD8_gate &
    Biotin_neg &
    ((Va2_neg & Vb5_neg) | Va2_xor_Vb5)
]

# Comp3 myeloid: untreated myeloid/DC
Myeloid_untreated_cells <- colnames(obj)[
  is_untreated &
    meta$final_clusters %in% myeloid_clusters
]

# Return to RNA for plotting
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# 5) Cluster label positions (median centers) + plotting helper
# ------------------------------------------------------------------------------
get_cluster_centers <- function(objx, reduction = "umap", cluster_col = "final_clusters") {
  emb <- Embeddings(objx, reduction = reduction)
  data.frame(
    x = emb[,1],
    y = emb[,2],
    cluster = as.character(objx[[cluster_col]][,1]),
    row.names = rownames(emb)
  ) %>%
    group_by(cluster) %>%
    summarise(x = median(x), y = median(y), .groups = "drop")
}

make_two_group_umap <- function(objx,
                                cells_myeloid,
                                cells_t,
                                title,
                                out_pdf,
                                reduction = "umap",
                                cluster_col = "final_clusters",
                                pt_size = 0.25) {
  grp <- rep("Other", ncol(objx))
  names(grp) <- colnames(objx)
  grp[colnames(objx) %in% cells_myeloid] <- "Myeloid"
  grp[colnames(objx) %in% cells_t] <- "T/OT"
  objx$plot_group <- factor(grp, levels = c("Other", "Myeloid", "T/OT"))

  centers <- get_cluster_centers(objx, reduction = reduction, cluster_col = cluster_col)

  p <- DimPlot(
    objx,
    reduction = reduction,
    group.by = "plot_group",
    pt.size = pt_size,
    cols = c(
      "Other" = "grey75",
      "Myeloid" = "#5fa9e3",  # darker pale blue
      "T/OT" = "#ea6d6d"   # darker pale red
    )
  ) +
    ggtitle(title) +
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

  pdf(out_pdf, width = 7.5, height = 6.5)
  print(p)
  dev.off()

  message("Saved: ", normalizePath(out_pdf, mustWork = FALSE))
}

# ------------------------------------------------------------------------------
# 6) THREE plots
# ------------------------------------------------------------------------------

# Plot 1: Comp1B (biotin only on OT-I; myeloid unfiltered)
make_two_group_umap(
  objx = obj,
  cells_myeloid = Myeloid_treated_cells,
  cells_t = OTI_comp1B_cells,
  title = "Comp1B (plot): OT-I (treated; biotin+) vs Treated Myeloid/DC (no biotin filter)",
  out_pdf = file.path(outdir, "UMAP_Comp1B_OTI_biotinOnlyOnOTI_vs_TreatedMyeloid.pdf")
)

# Plot 2: Comp2
make_two_group_umap(
  objx = obj,
  cells_myeloid = Myeloid_treated_cells,
  cells_t = EndogCD8_treated_cells,
  title = "Comp2: Endog CD8 (treated; CD8 ADT or RNA) vs Treated Myeloid/DC",
  out_pdf = file.path(outdir, "UMAP_Comp2_EndogCD8Treated_CD8_ADTorRNA_vs_TreatedMyeloid.pdf")
)

# Plot 3: Comp3
make_two_group_umap(
  objx = obj,
  cells_myeloid = Myeloid_untreated_cells,
  cells_t = EndogCD8_untreated_cells,
  title = "Comp3: Endog CD8 (untreated; CD8 ADT or RNA) vs Untreated Myeloid/DC",
  out_pdf = file.path(outdir, "UMAP_Comp3_EndogCD8Untreated_CD8_ADTorRNA_vs_UntreatedMyeloid.pdf")
)

# ------------------------------------------------------------------------------
# 7) Sanity counts
# ------------------------------------------------------------------------------
cat("\n--- Sanity counts (v2 fixed) ---\n")
cat("CD8 gate true (ADT or RNA):      ", sum(CD8_gate, na.rm = TRUE), "\n")
cat("Comp1B OT-I (treated; biotin+):  ", length(OTI_comp1B_cells), "\n")
cat("Comp2 Endog CD8 treated:         ", length(EndogCD8_treated_cells), "\n")
cat("Comp3 Endog CD8 untreated:       ", length(EndogCD8_untreated_cells), "\n")
cat("Myeloid treated (no biotin):     ", length(Myeloid_treated_cells), "\n")
cat("Myeloid untreated:               ", length(Myeloid_untreated_cells), "\n")
cat("-------------------------------\n")
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!exists("obj")) {
  seurat_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
  obj <- readRDS(seurat_rds)
}
stopifnot(exists("obj"))

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"),
                    "07_CELLCHAT_FEATURES_FINAL_BLUE")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Cluster label positions (robust)
# ------------------------------------------------------------------------------
get_cluster_centers <- function(objx) {
  emb <- Embeddings(objx, "umap")
  data.frame(
    x = emb[,1],
    y = emb[,2],
    cluster = as.character(objx$final_clusters)
  ) %>%
    group_by(cluster) %>%
    summarise(x = median(x), y = median(y), .groups="drop")
}

# ------------------------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------------------------
make_two_group_umap <- function(objx,
                                cells_myeloid,
                                cells_t,
                                title,
                                outfile) {

  grp <- rep("Other", ncol(objx))
  names(grp) <- colnames(objx)

  grp[colnames(objx) %in% cells_myeloid] <- "Myeloid"
  grp[colnames(objx) %in% cells_t] <- "T/OT"

  objx$plot_group <- factor(
    grp,
    levels=c("Other","Myeloid","T/OT")
  )

  centers <- get_cluster_centers(objx)

  p <- DimPlot(
    objx,
    reduction="umap",
    group.by="plot_group",
    pt.size=0.25,
    cols=c(
      "Other" = "grey75",
      "Myeloid" = "#859ee1",   # <- YOUR BLUE
      "T/OT" = "#ea6d6d"    # same darker red as before
    )
  ) +
    ggtitle(title) +
    theme_classic(base_size=12) +
    theme(
      plot.title = element_text(face="bold", hjust=0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    ) +
    geom_text(
      data=centers,
      aes(x=x,y=y,label=cluster),
      inherit.aes=FALSE,
      size=4,
      fontface="bold"
    )

  pdf(outfile, width=7.5, height=6.5)
  print(p)
  dev.off()
}

# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

# Comp1B
make_two_group_umap(
  obj,
  Myeloid_treated_cells,
  OTI_comp1B_cells,
  "Comp1B: OT-I (biotin+) vs Treated Myeloid/DC",
  file.path(outdir,"UMAP_Comp1B.pdf")
)

# Comp2
make_two_group_umap(
  obj,
  Myeloid_treated_cells,
  EndogCD8_treated_cells,
  "Comp2: Endogenous CD8 treated vs Treated Myeloid/DC",
  file.path(outdir,"UMAP_Comp2.pdf")
)

# Comp3
make_two_group_umap(
  obj,
  Myeloid_untreated_cells,
  EndogCD8_untreated_cells,
  "Comp3: Endogenous CD8 untreated vs Untreated Myeloid/DC",
  file.path(outdir,"UMAP_Comp3.pdf")
)

