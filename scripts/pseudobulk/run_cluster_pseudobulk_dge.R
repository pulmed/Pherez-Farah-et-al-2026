#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/pseudobulk/run_cluster_pseudobulk_dge.R
# Original file: DGE CLUSTERS (AVOID PSEUDOREPLICATION).R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Run paired sample-level pseudobulk edgeR differential-expression analyses between clusters.
# Inputs: Canonical analysis-ready Seurat object with RNA counts, final_clusters, and hash.ID metadata.
# Outputs: Pairwise edgeR CSV/XLSX tables and volcano PDFs.
# Assay/layer input: RNA counts.
# Dependencies: Seurat, edgeR, dplyr, tidyr, openxlsx, ggplot2, EnhancedVolcano.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - The paired design uses sample as the replication unit to avoid pseudoreplication.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggplot2)
  library(EnhancedVolcano)
})

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
RDS_FILE <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
OUTDIR <- Sys.getenv("OUTPUT_DIR", unset = file.path("output", "16_DGE_clusters_RNA_pseudobulk"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

CLUSTER_COL <- "final_clusters"
SAMPLE_COL <- "hash.ID"
ASSAY <- "RNA"

MIN_CELLS_PER_SAMPLE_PER_CLUSTER <- 20
MIN_SAMPLES_WITH_BOTH <- 2

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

MAKE_VOLCANO <- TRUE

# ------------------------------------------------------------------------------
# LOAD + CHECKS
# ------------------------------------------------------------------------------
if (!file.exists(RDS_FILE)) stop("Cannot find: ", RDS_FILE)
obj <- readRDS(RDS_FILE)

stopifnot(CLUSTER_COL %in% colnames(obj[[]]))
stopifnot(SAMPLE_COL  %in% colnames(obj[[]]))
stopifnot(ASSAY %in% names(obj@assays))

DefaultAssay(obj) <- ASSAY

counts <- GetAssayData(obj, assay = ASSAY, slot = "counts")

meta <- obj@meta.data %>%
  mutate(
    cell = rownames(.),
    cluster = as.character(.data[[CLUSTER_COL]]),
    sample = as.character(.data[[SAMPLE_COL]])
  ) %>%
  filter(!is.na(cluster), !is.na(sample))

clusters <- sort(unique(meta$cluster))

safe_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

# ------------------------------------------------------------------------------
# PSEUDOBULK BUILDER
# ------------------------------------------------------------------------------
make_pb_two_clusters <- function(cl1, cl2, meta, counts,
                                 min_cells = 20, min_samples = 2) {

  g1 <- paste0("cl_", cl1)
  g2 <- paste0("cl_", cl2)

  meta2 <- meta %>%
    filter(cluster %in% c(cl1, cl2)) %>%
    mutate(group = ifelse(cluster == cl1, g1, g2))

  tab <- meta2 %>%
    count(sample, group, name = "n_cells") %>%
    tidyr::pivot_wider(names_from = group, values_from = n_cells, values_fill = 0)

  # ensure columns exist even if one group absent
  if (!g1 %in% colnames(tab)) tab[[g1]] <- 0
  if (!g2 %in% colnames(tab)) tab[[g2]] <- 0

  keep_samples <- tab %>%
    filter(.data[[g1]] >= min_cells, .data[[g2]] >= min_cells) %>%
    pull(sample)

  if (length(keep_samples) < min_samples) {
    return(list(ok = FALSE, reason = "not_enough_samples_with_both_clusters",
                n_keep_samples = length(keep_samples)))
  }

  meta_keep <- meta2 %>% filter(sample %in% keep_samples)
  cells_keep <- meta_keep$cell
  sub_counts <- counts[, cells_keep, drop = FALSE]

  grp_fac <- interaction(meta_keep$sample, meta_keep$group, drop = TRUE, sep = "__")
  pb_counts <- t(rowsum(t(as.matrix(sub_counts)), group = grp_fac))

  split_mat <- do.call(rbind, strsplit(colnames(pb_counts), "__", fixed = TRUE))
  pheno <- data.frame(
    sample = factor(split_mat[, 1]),
    group = factor(split_mat[, 2], levels = c(g2, g1))  # baseline g2, test g1
  )
  rownames(pheno) <- colnames(pb_counts)

  chk <- table(pheno$sample, pheno$group)
  if (!all(chk[, g1] > 0 & chk[, g2] > 0)) {
    return(list(ok = FALSE, reason = "paired_design_not_satisfied_after_aggregation"))
  }

  list(ok = TRUE, pb_counts = pb_counts, pheno = pheno, tab_cells = tab, g1 = g1, g2 = g2)
}

# ------------------------------------------------------------------------------
# edgeR RUNNER (FIXED COEF)
# ------------------------------------------------------------------------------
run_edger_paired <- function(pb_counts, pheno) {
  dge <- DGEList(counts = pb_counts)

  design <- model.matrix(~ sample + group, data = pheno)
  # coefficient name for group second level vs first level
  coef_name <- paste0("group", levels(pheno$group)[2])

  if (!coef_name %in% colnames(design)) {
    stop("Coefficient not found in design matrix: ", coef_name,
         "\nDesign columns are: ", paste(colnames(design), collapse = ", "))
  }

  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = coef_name)

  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)
  list(tt = tt, design = design, coef_name = coef_name)
}

# ------------------------------------------------------------------------------
# MAIN LOOP (skip if already done)
# ------------------------------------------------------------------------------
for (i in 1:(length(clusters) - 1)) {
  cl1 <- clusters[i]
  message("[INFO] Comparing Cluster ", cl1)

  for (j in (i + 1):length(clusters)) {
    cl2 <- clusters[j]
    message("    -> Against Cluster ", cl2)

    base <- file.path(
      OUTDIR,
      paste0("Cluster", safe_filename(cl1), "_vs_Cluster", safe_filename(cl2))
    )

    out_xlsx <- paste0(base, "_edgeR.xlsx")
    out_pdf <- paste0(base, "_volcano.pdf")

    # DO NOT RUN TWICE
    if (file.exists(out_xlsx)) {
      message("    -> Already exists, skipping: ", basename(out_xlsx))
      next
    }

    pb <- make_pb_two_clusters(
      cl1 = cl1, cl2 = cl2,
      meta = meta, counts = counts,
      min_cells = MIN_CELLS_PER_SAMPLE_PER_CLUSTER,
      min_samples = MIN_SAMPLES_WITH_BOTH
    )

    if (!isTRUE(pb$ok)) {
      message("    [WARN] Skipping ", cl1, " vs ", cl2, " : ", pb$reason,
              " (n_keep_samples=", pb$n_keep_samples, ")")
      next
    }

    # Run edgeR
    res <- run_edger_paired(pb$pb_counts, pb$pheno)
    tt <- res$tt

    # logFC > 0 means higher in cl1 (g1) than cl2 (g2)
    tt$comparison <- paste0("Cluster ", cl1, " vs Cluster ", cl2)
    tt$direction <- ifelse(tt$logFC > 0, paste0("higher_in_", cl1), paste0("higher_in_", cl2))
    tt$sig_FDR <- tt$FDR <= FDR_CUTOFF
    tt$sig_FDR_logFC <- (tt$FDR <= FDR_CUTOFF) & (abs(tt$logFC) >= LOGFC_CUTOFF)

    tt_out <- tt %>%
      dplyr::select(comparison, gene, logFC, direction, logCPM, PValue, FDR, sig_FDR, sig_FDR_logFC)

    # Volcano
    if (MAKE_VOLCANO) {
      pdf(out_pdf, width = 7, height = 6)
      print(
        EnhancedVolcano(
          tt,
          lab = tt$gene,
          x = "logFC",
          y = "PValue",
          pCutoff = 0.05,
          FCcutoff = LOGFC_CUTOFF,
          title = paste0("Cluster ", cl1, " vs Cluster ", cl2),
          subtitle = paste0(
            "pseudobulk edgeR | design ~ sample + group\n",
            "coef: ", res$coef_name,
            " | samples used: ", length(unique(pb$pheno$sample)),
            " | min cells/sample/cluster: ", MIN_CELLS_PER_SAMPLE_PER_CLUSTER,
            "\n(logFC > 0 = higher in Cluster ", cl1, ")"
          ),
          pointSize = 1.6,
          labSize = 2.5
        )
      )
      dev.off()
    }

    # Save Excel (plus diagnostics)
    wb <- createWorkbook()

    addWorksheet(wb, "edgeR_results")
    writeData(wb, "edgeR_results", tt_out)

    addWorksheet(wb, "pseudobulk_cell_counts")
    writeData(wb, "pseudobulk_cell_counts", pb$tab_cells)

    addWorksheet(wb, "design_matrix")
    writeData(wb, "design_matrix", as.data.frame(res$design), rowNames = TRUE)

    addWorksheet(wb, "run_settings")
    writeData(wb, "run_settings", data.frame(
      rds_file = RDS_FILE,
      assay = ASSAY,
      cluster1 = cl1,
      cluster2 = cl2,
      group_baseline = pb$g2,
      group_tested = pb$g1,
      coef_tested = res$coef_name,
      min_cells_per_sample_per_cluster = MIN_CELLS_PER_SAMPLE_PER_CLUSTER,
      min_samples_with_both = MIN_SAMPLES_WITH_BOTH,
      fdr_cutoff = FDR_CUTOFF,
      logfc_cutoff = LOGFC_CUTOFF,
      stringsAsFactors = FALSE
    ))

    saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    # Optional CSV
    write.csv(tt_out, paste0(base, "_edgeR.csv"), row.names = FALSE)

    message("[OK] Saved pseudobulk edgeR for Cluster ", cl1, " vs ", cl2)
  }
}

message("Done. Output folder: ", OUTDIR)

