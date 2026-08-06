#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/visualization/cluster7_marker_heatmaps.R
# Original file: 20260217 HEATMAP 7 CELL AND SAMPLE.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Generate cluster 7 marker workbooks and cell/sample heatmaps.
# Inputs: Canonical analysis-ready Seurat object with RNA, SCT, final_clusters, and hash.ID metadata.
# Outputs: Cluster 7 top-marker XLSX files and cell/sample heatmap PDFs.
# Assay/layer input: RNA data and SCT data.
# Dependencies: Seurat, dplyr, openxlsx, pheatmap.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Runs both RNA and SCT assay variants.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(openxlsx)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
seu <- seurat_final_clean4

cluster_col <- "final_clusters"
target_cluster <- "7"

# Compare/visualize ONLY across these clusters, in this exact order:
cluster_order <- c("9","8","1","13","7")
clusters_compare <- cluster_order  # keep consistent

top_n <- 50
min_cells_per_sample_cluster <- 20
min_samples_per_cluster <- 2

candidate_sample_cols <- c("hashtag","HTO_maxID","sample","Sample","orig.ident","donor","Donor")

base_out_dir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "11_cluster7_top50_global_markers_compare_9_8_1_13_7")
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# HEATMAP APPEARANCE
# ------------------------------------------------------------------------------
hm_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100) # blue-white-red

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
detect_sample_col <- function(obj, candidates) {
  md <- obj@meta.data
  hit <- intersect(candidates, colnames(md))
  if (length(hit) == 0) stop("No sample/hashtag column detected in meta.data.")
  hit[1]
}

zscore_rows <- function(m) {
  mz <- t(scale(t(m)))
  mz[is.na(mz)] <- 0
  mz
}

avg_expr_one_assay <- function(obj, assay, features, group.by, slot = "data") {
  fn_formals <- names(formals(AverageExpression))
  if ("assays" %in% fn_formals) {
    AverageExpression(obj, assays = assay, features = features, group.by = group.by, slot = slot)[[assay]]
  } else {
    DefaultAssay(obj) <- assay
    AverageExpression(obj, features = features, group.by = group.by, slot = slot)[[1]]
  }
}

# handle Seurat numeric cluster naming: "7" vs "g7"
fix_cluster_cols <- function(mat, desired_ids) {
  if (all(desired_ids %in% colnames(mat))) return(desired_ids)
  g <- paste0("g", desired_ids)
  if (all(g %in% colnames(mat))) return(g)

  hit <- intersect(colnames(mat), c(desired_ids, g))
  if (length(hit) == 0) {
    stop("Could not match requested clusters to matrix columns.\nRequested: ",
         paste(desired_ids, collapse = ", "),
         "\nGot: ", paste(colnames(mat), collapse = ", "))
  }
  message("[WARN] Some clusters missing after filtering; using available: ", paste(hit, collapse = ", "))
  hit
}

cell_counts_table <- function(sample_vec, cluster_vec) {
  df <- as.data.frame(table(sample = sample_vec, cluster = cluster_vec), stringsAsFactors = FALSE)
  colnames(df)[3] <- "n_cells"
  df
}

ensure_rna_data_slot <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  # populate data slot if empty
  dat <- GetAssayData(obj, assay = "RNA", slot = "data")
  if (inherits(dat, "dgCMatrix") && length(dat@x) == 0) {
    message("-> Normalizing RNA to populate data slot")
    obj <- NormalizeData(obj, assay = "RNA")
  }
  obj
}

# ------------------------------------------------------------------------------
# CORE PIPELINE PER ASSAY
# ------------------------------------------------------------------------------
run_assay <- function(seu, assay = c("SCT","RNA")) {
  assay <- match.arg(assay)

  out_dir <- file.path(base_out_dir, assay)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_xlsx <- file.path(out_dir, paste0("cluster", target_cluster, "_top", top_n, "_", assay, ".xlsx"))
  hm_cell_pdf <- file.path(out_dir, paste0("HEATMAP_cell_level_", assay, ".pdf"))
  hm_samp_pdf <- file.path(out_dir, paste0("HEATMAP_sample_level_FILTERED_", assay, ".pdf"))

  # identities
  stopifnot(cluster_col %in% colnames(seu@meta.data))
  Idents(seu) <- cluster_col

  # sample column
  sample_col <- detect_sample_col(seu, candidate_sample_cols)

  # ------------------------------------------------------------------------------
  # 1) TOP 50 MARKERS: cluster 7 vs ALL clusters (global)
  # ------------------------------------------------------------------------------
  if (assay == "SCT") {
    stopifnot("SCT" %in% names(seu@assays))
    DefaultAssay(seu) <- "SCT"
    seu <- PrepSCTFindMarkers(seu)
  } else {
    stopifnot("RNA" %in% names(seu@assays))
    seu <- ensure_rna_data_slot(seu)
    DefaultAssay(seu) <- "RNA"
  }

  markers <- FindMarkers(
    seu,
    ident.1 = target_cluster,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  ) %>% rownames_to_column("gene")

  # filter mito / ribosomal (optional)
  markers_f <- markers %>%
    filter(
      !grepl("^mt-", gene, ignore.case = TRUE),
      !grepl("^(Rpl|Rps)", gene, ignore.case = TRUE)
    )

  lfc_col <- if ("avg_log2FC" %in% colnames(markers_f)) "avg_log2FC" else "avg_logFC"

  top_tbl <- markers_f %>%
    arrange(desc(.data[[lfc_col]])) %>%
    slice_head(n = top_n)

  genes_top <- top_tbl$gene

  # ------------------------------------------------------------------------------
  # 2) CELL-LEVEL COMPARISON: clusters 9,8,1,13,7 only
  # ------------------------------------------------------------------------------
  obj_sub <- subset(seu, idents = clusters_compare)
  DefaultAssay(obj_sub) <- assay

  avg_cell <- avg_expr_one_assay(
    obj_sub,
    assay = assay,
    features = genes_top,
    group.by = cluster_col,
    slot = "data"
  )

  keep_cols <- fix_cluster_cols(avg_cell, cluster_order)
  avg_cell <- avg_cell[, keep_cols, drop = FALSE]
  avg_cell_z <- zscore_rows(as.matrix(avg_cell))

  pdf(hm_cell_pdf, width = 8, height = 11)
  pheatmap(
    avg_cell_z,
    color = hm_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 12,
    fontface_row = "italic",
    border_color = NA,
    main = paste0("Top ", top_n, " cluster ", target_cluster, " markers (global) - ", assay)
  )
  dev.off()

  # ------------------------------------------------------------------------------
  # 3) SAMPLE-LEVEL (FILTERED) COMPARISON
  # ------------------------------------------------------------------------------
  md <- obj_sub@meta.data
  sample_vec <- as.character(md[[sample_col]])
  cluster_vec <- as.character(md[[cluster_col]])

  cc <- cell_counts_table(sample_vec, cluster_vec)
  valid_pairs <- cc[cc$n_cells >= min_cells_per_sample_cluster, , drop = FALSE]

  key_all <- paste(sample_vec, cluster_vec, sep = "__")
  key_valid <- paste(valid_pairs$sample, valid_pairs$cluster, sep = "__")
  cells_valid <- rownames(md)[key_all %in% key_valid]

  obj_pb <- subset(obj_sub, cells = cells_valid)
  DefaultAssay(obj_pb) <- assay

  # require >= min_samples_per_cluster per cluster
  md2 <- obj_pb@meta.data
  pairs <- unique(data.frame(
    sample = as.character(md2[[sample_col]]),
    cluster = as.character(md2[[cluster_col]]),
    stringsAsFactors = FALSE
  ))
  ns <- as.data.frame(table(pairs$cluster), stringsAsFactors = FALSE)
  colnames(ns) <- c("cluster","n_samples")
  valid_clusters <- ns$cluster[ns$n_samples >= min_samples_per_cluster]

  obj_pb <- subset(obj_pb, idents = intersect(valid_clusters, clusters_compare))
  DefaultAssay(obj_pb) <- assay

  # average per (sample, cluster)
  avg_sc <- avg_expr_one_assay(
    obj_pb,
    assay = assay,
    features = genes_top,
    group.by = c(sample_col, cluster_col),
    slot = "data"
  )

  # parse columns like SAMPLE_CLUSTER (cluster may be g7)
  cn <- colnames(avg_sc)
  cl_parsed <- sub("^.*_([^_]*)$", "\\1", cn)
  samp_parsed <- sub("^(.*)_[^_]*$", "\\1", cn)
  col_info <- data.frame(col = cn, sample = samp_parsed, cluster = cl_parsed, stringsAsFactors = FALSE)

  avg_long <- as.data.frame(avg_sc) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "col", values_to = "expr") %>%
    left_join(col_info, by = "col")

  # mean across samples within each cluster
  avg_samp <- avg_long %>%
    group_by(gene, cluster) %>%
    summarise(mean_expr = mean(expr), .groups = "drop") %>%
    pivot_wider(names_from = cluster, values_from = mean_expr) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  keep_cols2 <- fix_cluster_cols(avg_samp, cluster_order)
  avg_samp <- avg_samp[, keep_cols2, drop = FALSE]
  avg_samp_z <- zscore_rows(avg_samp)

  pdf(hm_samp_pdf, width = 8, height = 11)
  pheatmap(
    avg_samp_z,
    color = hm_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 12,
    fontface_row = "italic",
    border_color = NA,
    main = paste0("Sample-level mean (filtered) - ", assay)
  )
  dev.off()

  # ------------------------------------------------------------------------------
  # 4) EXPORT EXCEL
  # ------------------------------------------------------------------------------
  wb <- createWorkbook()

  addWorksheet(wb, paste0("Top", top_n, "_markers_cluster", target_cluster, "_global"))
  writeData(wb, 1, top_tbl)

  addWorksheet(wb, "CellLevel_clusterAvg")
  writeData(wb, "CellLevel_clusterAvg",
            cbind(gene = rownames(avg_cell), as.data.frame(avg_cell)))

  addWorksheet(wb, "CellLevel_rowZ")
  writeData(wb, "CellLevel_rowZ",
            cbind(gene = rownames(avg_cell_z), as.data.frame(avg_cell_z)))

  addWorksheet(wb, "SampleLevel_clusterMean")
  writeData(wb, "SampleLevel_clusterMean",
            cbind(gene = rownames(avg_samp), as.data.frame(avg_samp)))

  addWorksheet(wb, "SampleLevel_rowZ")
  writeData(wb, "SampleLevel_rowZ",
            cbind(gene = rownames(avg_samp_z), as.data.frame(avg_samp_z)))

  addWorksheet(wb, "Filtering_info")
  writeData(wb, "Filtering_info", data.frame(
    parameter = c("assay_used",
                  "cluster_col",
                  "sample_col_used",
                  "target_cluster",
                  "clusters_compared_order",
                  "top_n",
                  "min_cells_per_sample_cluster",
                  "min_samples_per_cluster"),
    value = c(assay,
              cluster_col,
              sample_col,
              target_cluster,
              paste(cluster_order, collapse = ","),
              top_n,
              min_cells_per_sample_cluster,
              min_samples_per_cluster)
  ))

  saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  message("[OK] Finished ", assay)
  message("   Folder:  ", out_dir)
  message("   Excel:   ", out_xlsx)
  message("   Cell HM: ", hm_cell_pdf)
  message("   Samp HM: ", hm_samp_pdf)
}

# ------------------------------------------------------------------------------
# RUN BOTH ASSAYS
# ------------------------------------------------------------------------------
run_assay(seu, "SCT")
run_assay(seu, "RNA")

message("\n[OK] All outputs saved under:\n", base_out_dir)

