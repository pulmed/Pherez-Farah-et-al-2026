#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/cellchat/run_biotin_comparisons.R
# Original file: CellChat v2.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Prepare biotin-aware groupings and run CellChat communication analyses.
# Inputs: Canonical analysis-ready Seurat object with SCT, ADT, final_clusters, and condition metadata.
# Outputs: CellChat objects, communication tables, UMAP group plots, and chord/interaction plots.
# Assay/layer input: SCT data for CellChat expression; ADT data for Biotin and CD8 gates.
# Dependencies: Seurat, CellChat, dplyr, openxlsx, ggplot2, tidyverse, svglite, scales.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Set SMOKE_TEST=1 to validate inputs without running full CellChat analysis.
# - CD8 RNA helper flags are calculated, but the final CD8_pos gate uses ADT CD8a-TotalSeqC.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

library(Seurat)
library(CellChat)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(svglite)
library(scales)

options(stringsAsFactors = FALSE)

# Input / output
seurat_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
outdir_root <- Sys.getenv("OUTPUT_DIR", unset = file.path("output", "06_cellchat_biotin_normalized"))
dir.create(outdir_root, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load Seurat object & basic checks
# ------------------------------------------------------------------------------
seurat_object <- readRDS(seurat_rds)

stopifnot("final_clusters" %in% colnames(seurat_object@meta.data))
stopifnot("condition"      %in% colnames(seurat_object@meta.data))

if (identical(Sys.getenv("SMOKE_TEST"), "1")) {
  message("SMOKE_TEST=1: input object and required CellChat metadata validated; skipping full CellChat analysis.")
  quit(save = "no", status = 0)
}

Idents(seurat_object) <- seurat_object$final_clusters
meta <- seurat_object@meta.data

# ------------------------------------------------------------------------------
# 2. Compute biotin_status (cluster-specific thresholds)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 2A. Normalize ADT (CLR) & helper to extract biotin ADT
# ------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "ADT"
seurat_object <- NormalizeData(seurat_object, normalization.method = "CLR", margin = 2)

get_biotin <- function(obj, cells) {
  GetAssayData(obj, assay = "ADT", slot = "data")["Biotin-TotalSeqC", cells, drop = TRUE]
}

clusters <- sort(unique(seurat_object$final_clusters))
thresholds <- list()
untreated_counts <- integer(length(clusters))
names(untreated_counts) <- as.character(clusters)

# ------------------------------------------------------------------------------
# 2B. Per-cluster thresholds from UNTREATED cells (Q3)
# ------------------------------------------------------------------------------
for (cl in clusters) {
  cells_u <- WhichCells(
    seurat_object,
    expression = final_clusters == cl & condition == "untreated"
  )

  untreated_counts[as.character(cl)] <- length(cells_u)

  if (length(cells_u) >= 5) {
    vals <- get_biotin(seurat_object, cells_u)
    thresholds[[as.character(cl)]] <- as.numeric(quantile(vals, 0.75, na.rm = TRUE))
  } else {
    thresholds[[as.character(cl)]] <- NA_real_
  }
}

thresholds <- unlist(thresholds)

# Optional: inspect thresholds
biotin_thr_df <- data.frame(
  cluster = names(untreated_counts),
  n_untreated = as.integer(untreated_counts),
  Q3 = thresholds[match(names(untreated_counts), names(thresholds))]
)
biotin_thr_df

# ------------------------------------------------------------------------------
# 2C. Assign biotin_status to ALL cells
# ------------------------------------------------------------------------------
biotin_status <- rep(NA_character_, ncol(seurat_object))
names(biotin_status) <- colnames(seurat_object)

for (cl in names(thresholds)) {
  thr <- thresholds[[cl]]
  if (is.na(thr)) next

  cells_cl <- WhichCells(seurat_object, expression = final_clusters == as.numeric(cl))
  if (length(cells_cl) == 0) next

  vals <- get_biotin(seurat_object, cells_cl)
  biotin_status[names(vals)] <- ifelse(vals > thr, "positive", "negative")
}

seurat_object$biotin_status <- factor(biotin_status, levels = c("negative", "positive"))

# Return to SCT for expression & CellChat
DefaultAssay(seurat_object) <- "RNA"

# ------------------------------------------------------------------------------
# 3. Major group mapping (your 7 groups)
# ------------------------------------------------------------------------------

cluster_groups <- list(
  `1` = c(0,2,3,4,5,14,17),  # Main myeloid
  `2` = c(1,7,8,9,13),             # Dendritic cells
  `3` = c(6,10,11,15),          # T cell cluster
  `4` = 12,                     # Tumor cluster
  `5` = 16,                     # pDCs
  `6` = 18,                     # CAFs
  `7` = 19                      # TECs
)

major_group <- rep(NA_character_, ncol(seurat_object))
names(major_group) <- colnames(seurat_object)

fc <- seurat_object$final_clusters

for (grp in names(cluster_groups)) {
  cl_vec <- cluster_groups[[grp]]
  cells_grp <- colnames(seurat_object)[fc %in% cl_vec]
  major_group[cells_grp] <- switch(
    grp,
    `1` = "Myeloid_main",
    `2` = "Dendritic_cells",
    `3` = "T_cells",
    `4` = "Tumor",
    `5` = "pDCs",
    `6` = "CAFs",
    `7` = "TECs"
  )
}

seurat_object$major_group <- factor(
  major_group,
  levels = c("Myeloid_main","Dendritic_cells","T_cells","Tumor","pDCs","CAFs","TECs")
)

# ------------------------------------------------------------------------------
# 4. Define OT-I, endogenous CD8 T cells, and myeloid cells
# ------------------------------------------------------------------------------

adt <- GetAssayData(seurat_object, assay = "ADT", slot = "data")
rna <- GetAssayData(seurat_object, assay = "SCT", slot = "data")
DefaultAssay(seurat_object) <- "SCT"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
feature_pos <- function(mat, feature, threshold = 0) {
  if (!feature %in% rownames(mat)) {
    return(rep(FALSE, ncol(mat)))
  }
  vals <- mat[feature, , drop = TRUE]
  vals > threshold
}

# ADT positivity
Va2_pos <- feature_pos(adt, "Va2-TotalSeqC")
Vb5_pos <- feature_pos(adt, "Vb5-TotalSeqC")
CD8a_adt_pos <- feature_pos(adt, "CD8a-TotalSeqC")
Biotin_raw_pos <- feature_pos(adt, "Biotin-TotalSeqC")

# RNA positivity
Cd8a_rna_pos <- feature_pos(rna, "Cd8a")
Cd8b1_rna_pos <- feature_pos(rna, "Cd8b1")

CD8_pos <- CD8a_adt_pos

Biotin_pos <- seurat_object$biotin_status == "positive"
Biotin_neg <- seurat_object$biotin_status == "negative"

Va2_neg <- !Va2_pos
Vb5_neg <- !Vb5_pos
Biotin_raw_neg <- !Biotin_raw_pos

# XOR: Va2 or Vb5, but not both
Va2_xor_Vb5 <- (Va2_pos & !Vb5_pos) | (!Va2_pos & Vb5_pos)

# T-cell clusters now ONLY 6,10,11 (exclude 15)
tcell_clusters <- c(6,10,11)

# Final myeloid/DC clusters of interest
myeloid_clusters <- c(0,1,2,3,4,5,7,8,9,13,14,16,17)

# Convenience flags
is_treated <- meta$condition == "treated"
is_untreated <- meta$condition == "untreated"

# ------------------------------------------------------------------------------
# 4A. OT-I: treated, T-cell clusters, Va2+ & Vb5+, CD8a+, biotin positive
# ------------------------------------------------------------------------------
OTI_cells <- colnames(seurat_object)[
  is_treated &
  meta$final_clusters %in% tcell_clusters &
  Va2_pos & Vb5_pos &
  CD8_pos &
  Biotin_pos
]

# ------------------------------------------------------------------------------
# 4B. Endogenous CD8 T cells (treated): CD8a+, biotin negative
# ------------------------------------------------------------------------------
EndogCD8_treated_cells <- colnames(seurat_object)[
  is_treated &
  meta$final_clusters %in% tcell_clusters &
  CD8_pos &
  Biotin_neg
]

# ------------------------------------------------------------------------------
# 4C. Endogenous CD8 T cells (untreated)
# ------------------------------------------------------------------------------
EndogCD8_untreated_cells <- colnames(seurat_object)[
  is_untreated &
  meta$final_clusters %in% tcell_clusters &
  CD8_pos &
  Biotin_neg
]

# ------------------------------------------------------------------------------
# 4D. Myeloid/DC cells for different uses
# ------------------------------------------------------------------------------
# A-type (no biotin filter), treated, hashtags 1-4
Myeloid_treated_A_cells <- colnames(seurat_object)[
  is_treated &
  meta$final_clusters %in% myeloid_clusters
]

# B-type (biotin threshold; biotin_status == "positive"), treated, hashtags 1-4
Myeloid_treated_B_cells <- colnames(seurat_object)[
  is_treated &
  meta$final_clusters %in% myeloid_clusters &
  Biotin_pos
]

# Endogenous treated comparison: same as A
Myeloid_treated_for_endog_cells <- Myeloid_treated_A_cells

# Untreated myeloid/DC: hashtags 5-8
Myeloid_untreated_cells <- colnames(seurat_object)[
  is_untreated &
  meta$final_clusters %in% myeloid_clusters
]

# Sanity checks
cat("OT-I cells:                      ", length(OTI_cells), "\n")
cat("Endog CD8 treated cells:         ", length(EndogCD8_treated_cells), "\n")
cat("Endog CD8 untreated cells:       ", length(EndogCD8_untreated_cells), "\n")
cat("Myeloid treated (A, no biotin):  ", length(Myeloid_treated_A_cells), "\n")
cat("Myeloid treated (B, biotin+):    ", length(Myeloid_treated_B_cells), "\n")
cat("Myeloid untreated:               ", length(Myeloid_untreated_cells), "\n")

# ------------------------------------------------------------------------------
# 5. UMAP sanity plots: major_group & comm_groups
# ------------------------------------------------------------------------------

# 5A. Global UMAP with major_group
pdf(file.path(outdir_root, "UMAP_major_group.pdf"), width = 7, height = 6)
print(
  DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = "major_group",
    label = TRUE,
    pt.size = 0.3
  )
)
dev.off()

# 5B. UMAP for communication groups (OTI / Endog / Myeloid / Other)

seurat_object$comm_groups <- "Other"

seurat_object$comm_groups[colnames(seurat_object) %in% OTI_cells] <- "OTI_Tcell"
seurat_object$comm_groups[colnames(seurat_object) %in% EndogCD8_treated_cells] <- "EndogCD8_T_treated"
seurat_object$comm_groups[colnames(seurat_object) %in% EndogCD8_untreated_cells] <- "EndogCD8_T_untreated"
seurat_object$comm_groups[colnames(seurat_object) %in% Myeloid_treated_A_cells] <- "Myeloid_treated"
seurat_object$comm_groups[colnames(seurat_object) %in% Myeloid_untreated_cells] <- "Myeloid_untreated"

seurat_object$comm_groups <- factor(
  seurat_object$comm_groups,
  levels = c("OTI_Tcell",
             "EndogCD8_T_treated",
             "EndogCD8_T_untreated",
             "Myeloid_treated",
             "Myeloid_untreated",
             "Other")
)

pdf(file.path(outdir_root, "UMAP_comm_groups_allcells.pdf"), width = 7, height = 6)
print(
  DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = "comm_groups",
    label = TRUE,
    pt.size = 0.3
  )
)
dev.off()

# Only cells that participate in CellChat comparisons
cells_used <- unique(c(
  OTI_cells,
  EndogCD8_treated_cells,
  EndogCD8_untreated_cells,
  Myeloid_treated_A_cells,
  Myeloid_treated_B_cells,
  Myeloid_untreated_cells
))

pdf(file.path(outdir_root, "UMAP_comm_groups_cells_used.pdf"), width = 7, height = 6)
print(
  DimPlot(
    subset(seurat_object, cells = cells_used),
    reduction = "umap",
    group.by = "comm_groups",
    label = TRUE,
    pt.size = 0.5
  )
)
dev.off()

# ------------------------------------------------------------------------------
# 6. Helper: run CellChat & export LR tables (csv + xlsx)
# ------------------------------------------------------------------------------

run_cellchat_and_export <- function(
  seurat_sub,
  group_column = "celltype_communication",
  t_label,
  species = "mouse",
  comparison_label,
  outdir,
  min_cells = 10,
  sample_col = "hash.ID",   # e.g. "orig.ident" or "hashtag" or "sample"
  min_samples_per_grp = 2
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Use RNA assay for CellChat
  DefaultAssay(seurat_sub) <- "SCT"
  data.input <- GetAssayData(seurat_sub, assay = "SCT", slot = "data")
  meta.sub <- seurat_sub@meta.data

  # ------------------------------------------------------------------------------
  # Auto-detect sample column if not provided
  # ------------------------------------------------------------------------------
  if (is.null(sample_col)) {
    candidates <- c("orig.ident", "hashtag", "sample", "Sample", "donor", "Donor", "replicate", "Replicate")
    sample_col <- candidates[candidates %in% colnames(meta.sub)][1]
    if (is.na(sample_col) || is.null(sample_col)) {
      stop("No sample column found. Please pass sample_col = 'orig.ident' (or your sample id column).")
    }
  }
  stopifnot(sample_col %in% colnames(meta.sub))
  stopifnot(group_column %in% colnames(meta.sub))

  # ------------------------------------------------------------------------------
  # Enforce group-level thresholds:
  #   - at least min_cells cells
  #   - at least min_samples_per_grp distinct samples
  # ------------------------------------------------------------------------------
  qc <- meta.sub %>%
    dplyr::select(all_of(group_column), all_of(sample_col)) %>%
    dplyr::mutate(
      grp = .data[[group_column]],
      samp = .data[[sample_col]]
    ) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(samp),
      .groups = "drop"
    ) %>%
    dplyr::arrange(grp)

  # save QC so you can show Alfredo quickly
  qc_xlsx <- file.path(outdir, paste0("QC_group_counts_", comparison_label, ".xlsx"))
  wb_qc <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_qc, "QC")
  openxlsx::writeData(wb_qc, "QC", qc)
  openxlsx::saveWorkbook(wb_qc, qc_xlsx, overwrite = TRUE)

  keep_grps <- qc %>%
    dplyr::filter(n_cells >= min_cells, n_samples >= min_samples_per_grp) %>%
    dplyr::pull(grp)

  keep_cells <- rownames(meta.sub)[meta.sub[[group_column]] %in% keep_grps]

  # Subset to kept cells before running CellChat
  data.input <- data.input[, keep_cells, drop = FALSE]
  meta.sub <- meta.sub[keep_cells, , drop = FALSE]
  meta.sub[[group_column]] <- droplevels(factor(meta.sub[[group_column]]))

  # ------------------------------------------------------------------------------
  # Create CellChat object
  # ------------------------------------------------------------------------------
  cellchat <- createCellChat(object = data.input, meta = meta.sub, group.by = group_column)

  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else {
    CellChatDB <- CellChatDB.human
  }
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"))
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)

  # Still keep the CellChat filter (works on final communication table)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  df_comm <- subsetCommunication(cellchat)

  # Identify myeloid groups (Myeloid_<clusterid>)
  myeloid_groups <- grep("^Myeloid_", unique(c(df_comm$source, df_comm$target)), value = TRUE)

  # Filter to ONLY T <-> Myeloid interactions (both directions)
  df_T_myeloid <- df_comm %>%
    dplyr::filter(
      (source == t_label & target %in% myeloid_groups) |
      (target == t_label & source %in% myeloid_groups)
    )

  # Excel with 1 sheet per myeloid cluster + overview
  xlsx_file <- file.path(outdir, paste0("CellChat_LR_", comparison_label, "_T_vs_Myeloid.xlsx"))
  wb <- createWorkbook()

  addWorksheet(wb, "All_T_vs_Myeloid")
  writeData(wb, sheet = "All_T_vs_Myeloid", x = df_T_myeloid)

  for (mg in myeloid_groups) {
    df_mg <- df_T_myeloid %>%
      dplyr::filter(
        (source == t_label & target == mg) |
          (source == mg      & target == t_label)
      )
    if (nrow(df_mg) == 0) next

    sheet_name <- gsub("Myeloid_", "cl_", mg)
    sheet_name <- substr(sheet_name, 1, 30)

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = df_mg)
  }

  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  message("Saved Excel (T vs Myeloid only): ", xlsx_file)
  message("Saved QC: ", qc_xlsx)

  return(cellchat)
}

# ------------------------------------------------------------------------------
# 7. COMPARISON 1A:
#    OT-I vs Treated Myeloid (hashtags 1-4, no biotin filter)
# ------------------------------------------------------------------------------

cells_comp1A <- union(OTI_cells, Myeloid_treated_A_cells)
obj_comp1A <- subset(seurat_object, cells = cells_comp1A)

obj_comp1A$celltype_communication <- NA_character_
obj_comp1A$celltype_communication[colnames(obj_comp1A) %in% OTI_cells] <- "OTI_Tcell"

myeloid_idx_1A <- colnames(obj_comp1A) %in% Myeloid_treated_A_cells
obj_comp1A$celltype_communication[myeloid_idx_1A] <-
  paste0("Myeloid_", obj_comp1A$final_clusters[myeloid_idx_1A])

obj_comp1A$celltype_communication <- factor(obj_comp1A$celltype_communication)
Idents(obj_comp1A) <- "celltype_communication"

outdir_comp1A <- file.path(outdir_root, "Comp1A_OTI_vs_TreatedMyeloid_noBiotin")
cellchat_comp1A <- run_cellchat_and_export(
  seurat_sub = obj_comp1A,
  group_column = "celltype_communication",
  t_label = "OTI_Tcell",
  comparison_label = "OTI_vs_TreatedMyeloid_noBiotin",
  outdir = outdir_comp1A
)

# ------------------------------------------------------------------------------
# 8. COMPARISON 1B:
#    OT-I vs Treated Myeloid_Biotin+ (hashtags 1-4, myeloid biotin+)
# ------------------------------------------------------------------------------

cells_comp1B <- union(OTI_cells, Myeloid_treated_B_cells)
obj_comp1B <- subset(seurat_object, cells = cells_comp1B)

obj_comp1B$celltype_communication <- NA_character_
obj_comp1B$celltype_communication[colnames(obj_comp1B) %in% OTI_cells] <- "OTI_Tcell"

myeloid_idx_1B <- colnames(obj_comp1B) %in% Myeloid_treated_B_cells
obj_comp1B$celltype_communication[myeloid_idx_1B] <-
  paste0("Myeloid_", obj_comp1B$final_clusters[myeloid_idx_1B])

obj_comp1B$celltype_communication <- factor(obj_comp1B$celltype_communication)
Idents(obj_comp1B) <- "celltype_communication"

outdir_comp1B <- file.path(outdir_root, "Comp1B_OTI_vs_TreatedMyeloid_BiotinPos")
cellchat_comp1B <- run_cellchat_and_export(
  seurat_sub = obj_comp1B,
  group_column = "celltype_communication",
  t_label = "OTI_Tcell",
  comparison_label = "OTI_vs_TreatedMyeloid_BiotinPos",
  outdir = outdir_comp1B
)

# ------------------------------------------------------------------------------
# 9. COMPARISON 2:
#    Endogenous CD8a (treated) vs Treated Myeloid (A definition)
# ------------------------------------------------------------------------------

cells_comp2 <- union(EndogCD8_treated_cells, Myeloid_treated_for_endog_cells)
obj_comp2 <- subset(seurat_object, cells = cells_comp2)

obj_comp2$celltype_communication <- NA_character_
obj_comp2$celltype_communication[colnames(obj_comp2) %in% EndogCD8_treated_cells] <-
  "EndogCD8_T_treated"

myeloid_idx_2 <- colnames(obj_comp2) %in% Myeloid_treated_for_endog_cells
obj_comp2$celltype_communication[myeloid_idx_2] <-
  paste0("Myeloid_", obj_comp2$final_clusters[myeloid_idx_2])

obj_comp2$celltype_communication <- factor(obj_comp2$celltype_communication)
Idents(obj_comp2) <- "celltype_communication"

outdir_comp2 <- file.path(outdir_root, "Comp2_EndogCD8_Treated_vs_TreatedMyeloid")
cellchat_comp2 <- run_cellchat_and_export(
  seurat_sub = obj_comp2,
  group_column = "celltype_communication",
  t_label = "EndogCD8_T_treated",
  comparison_label = "EndogCD8_Treated_vs_TreatedMyeloid",
  outdir = outdir_comp2
)

# ------------------------------------------------------------------------------
# 10. COMPARISON 3:
#     Endogenous CD8a (untreated) vs Untreated Myeloid
# ------------------------------------------------------------------------------

cells_comp3 <- union(EndogCD8_untreated_cells, Myeloid_untreated_cells)
obj_comp3 <- subset(seurat_object, cells = cells_comp3)

obj_comp3$celltype_communication <- NA_character_
obj_comp3$celltype_communication[colnames(obj_comp3) %in% EndogCD8_untreated_cells] <-
  "EndogCD8_T_untreated"

myeloid_idx_3 <- colnames(obj_comp3) %in% Myeloid_untreated_cells
obj_comp3$celltype_communication[myeloid_idx_3] <-
  paste0("Myeloid_", obj_comp3$final_clusters[myeloid_idx_3])

obj_comp3$celltype_communication <- factor(obj_comp3$celltype_communication)
Idents(obj_comp3) <- "celltype_communication"

outdir_comp3 <- file.path(outdir_root, "Comp3_EndogCD8_Untreated_vs_UntreatedMyeloid")
cellchat_comp3 <- run_cellchat_and_export(
  seurat_sub = obj_comp3,
  group_column = "celltype_communication",
  t_label = "EndogCD8_T_untreated",
  comparison_label = "EndogCD8_Untreated_vs_UntreatedMyeloid",
  outdir = outdir_comp3
)

# ------------------------------------------------------------------------------
# 11. Chord gene plots: per myeloid cluster, per comparison
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Global style knobs
# ------------------------------------------------------------------------------
LAB_CEX <- 0.9
PDF_SIZE <- 14
TRANSPARENCY <- 0.15

# ------------------------------------------------------------------------------
# Fixed OT-I color
# ------------------------------------------------------------------------------
COL_OTI <- "#E09E9A"
COL_MYELOID <- "#7792CA"

# ------------------------------------------------------------------------------
# Main plotting function: all myeloid clusters
# ------------------------------------------------------------------------------
make_chord_gene_per_myeloid <- function(cellchat_obj,
                                        sender_name,
                                        comparison_outdir,
                                        title_prefix = "") {

  df_comm <- subsetCommunication(cellchat_obj) %>%
    dplyr::filter(source != target) %>%
    dplyr::filter(!is.na(ligand), !is.na(receptor)) %>%
    dplyr::filter(ligand != receptor)

  myeloid_groups <- grep("^Myeloid_", unique(c(df_comm$source, df_comm$target)), value = TRUE)

  output_dir <- file.path(comparison_outdir, "Chords_PerCluster_2nodes_gene")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (cluster in sort(myeloid_groups)) {
    message("[", sender_name, "] chord plot for: ", cluster)

    groups <- c(sender_name, cluster)

    df_subset <- df_comm %>%
      dplyr::filter(
        (source == sender_name & target == cluster) |
        (source == cluster & target == sender_name)
      )

    if (nrow(df_subset) == 0) {
      message("Skipping ", cluster, ": no valid interactions.")
      next
    }

    original_comm <- cellchat_obj@net$communication
    cellchat_obj@net$communication <- df_subset

    color_map <- c(
      setNames(COL_OTI, sender_name),
      setNames(COL_MYELOID, cluster)
    )

    base_name <- paste0("ChordGene_", title_prefix, sender_name, "_vs_", cluster)
    pdf_filename <- file.path(output_dir, paste0(base_name, ".pdf"))
    svg_filename <- file.path(output_dir, paste0(base_name, ".svg"))

    pdf(pdf_filename, width = PDF_SIZE, height = PDF_SIZE, useDingbats = FALSE)
    tryCatch({
      netVisual_chord_gene(
        object = cellchat_obj,
        sources.use = groups,
        targets.use = groups,
        net = df_subset,
        color.use = color_map[groups],
        lab.cex = LAB_CEX,
        small.gap = 4,
        big.gap = 8,
        transparency = TRANSPARENCY,
        title.name = paste0(title_prefix, sender_name, " <-> ", cluster),
        show.legend = TRUE
      )
    }, error = function(e) {
      message("Skipping ", cluster, " due to error: ", e$message)
    })
    dev.off()

    svglite(svg_filename, width = PDF_SIZE, height = PDF_SIZE)
    tryCatch({
      netVisual_chord_gene(
        object = cellchat_obj,
        sources.use = groups,
        targets.use = groups,
        net = df_subset,
        color.use = color_map[groups],
        lab.cex = LAB_CEX,
        small.gap = 4,
        big.gap = 8,
        transparency = TRANSPARENCY,
        title.name = paste0(title_prefix, sender_name, " <-> ", cluster),
        show.legend = TRUE
      )
    }, error = function(e) {
      message("Skipping ", cluster, " due to error: ", e$message)
    })
    dev.off()

    cellchat_obj@net$communication <- original_comm
  }
}

# 11A. Chords for Comp1A
make_chord_gene_per_myeloid(
  cellchat_obj = cellchat_comp1A,
  sender_name = "OTI_Tcell",
  comparison_outdir = outdir_comp1A,
  title_prefix = "Comp1A_OTI_"
)

# 11B. Chords for Comp1B
make_chord_gene_per_myeloid(
  cellchat_obj = cellchat_comp1B,
  sender_name = "OTI_Tcell",
  comparison_outdir = outdir_comp1B,
  title_prefix = "Comp1B_OTI_BiotinPos_"
)

# 11C. Chords for Comp2
make_chord_gene_per_myeloid(
  cellchat_obj = cellchat_comp2,
  sender_name = "EndogCD8_T_treated",
  comparison_outdir = outdir_comp2,
  title_prefix = "Comp2_EndogCD8_Treated_"
)

# 11D. Chords for Comp3
make_chord_gene_per_myeloid(
  cellchat_obj = cellchat_comp3,
  sender_name = "EndogCD8_T_untreated",
  comparison_outdir = outdir_comp3,
  title_prefix = "Comp3_EndogCD8_Untreated_"
)

# ------------------------------------------------------------------------------
# 12. UNIQUE INTERACTIONS exports
#     - File #1: unique interactions in Comp1A vs (Comp2 + Comp3)
#     - File #2: unique interactions in Comp1B vs (Comp2 + Comp3)
#     - File #3: union of (File #1 and File #2) with a column indicating:
#           only_1A_vs_2plus3 / only_1B_vs_2plus3 / both_1A_and_1B_vs_2plus3
#     - Sheets: batch_<cluster> for Monocyte/Macrophage + DC clusters
# ------------------------------------------------------------------------------
# 12A. Input excel paths
# ------------------------------------------------------------------------------
# These Excel files are produced earlier by run_cellchat_and_export()
# (make sure these filenames match what you actually wrote out)

comp1A_xlsx <- file.path(
  outdir_comp1A,
  "CellChat_LR_OTI_vs_TreatedMyeloid_noBiotin_T_vs_Myeloid.xlsx"
)

comp1B_xlsx <- file.path(
  outdir_comp1B,
  "CellChat_LR_OTI_vs_TreatedMyeloid_BiotinPos_T_vs_Myeloid.xlsx"
)

comp2_xlsx <- file.path(
  outdir_comp2,
  "CellChat_LR_EndogCD8_Treated_vs_TreatedMyeloid_T_vs_Myeloid.xlsx"
)

comp3_xlsx <- file.path(
  outdir_comp3,
  "CellChat_LR_EndogCD8_Untreated_vs_UntreatedMyeloid_T_vs_Myeloid.xlsx"
)

# ------------------------------------------------------------------------------
# 12B. Cluster sets to keep
# ------------------------------------------------------------------------------
mono_mac_clusters <- c(0,2,3,4,5,14,17)  # Myeloid_main
dc_clusters <- c(1,7,8,9,13)       # Dendritic_cells
clusters_keep <- sort(c(mono_mac_clusters, dc_clusters))

# ------------------------------------------------------------------------------
# 12C. Helpers
# ------------------------------------------------------------------------------
read_all_sheet <- function(xlsx_path, sheet = "All_T_vs_Myeloid") {
  if (!file.exists(xlsx_path)) {
    stop("File not found: ", xlsx_path)
  }
  df <- read.xlsx(xlsx_path, sheet = sheet)

  needed <- c("source","target","ligand","receptor")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(xlsx_path), ": ",
         paste(missing, collapse = ", "))
  }

  # CellChat version differences: ensure prob exists
  if (!("prob" %in% colnames(df))) df$prob <- NA_real_

  df
}

extract_myeloid_cluster <- function(source_vec, target_vec) {
  src_cl <- suppressWarnings(as.integer(str_match(source_vec, "^Myeloid_(\\d+)$")[,2]))
  tgt_cl <- suppressWarnings(as.integer(str_match(target_vec, "^Myeloid_(\\d+)$")[,2]))
  ifelse(!is.na(src_cl), src_cl, tgt_cl)
}

# Determine direction and role relative to the T population.
# In your exports, the T label is always in the source/target, so this is robust.
infer_direction_role <- function(source_vec, target_vec) {
  # If source is NOT Myeloid_* then it's the T label (OTI_Tcell or EndogCD8_*),
  # so call it T_to_Myeloid. Else Myeloid_to_T.
  source_is_myeloid <- grepl("^Myeloid_\\d+$", source_vec)
  direction <- ifelse(source_is_myeloid, "Myeloid_to_T", "T_to_Myeloid")
  role <- ifelse(direction == "T_to_Myeloid", "Outgoing_from_T", "Incoming_to_T")
  list(direction = direction, role = role)
}

prep_df <- function(df) {
  dr <- infer_direction_role(df$source, df$target)

  df %>%
    mutate(
      myeloid_cluster = extract_myeloid_cluster(source, target),
      Direction = dr$direction,
      Role = dr$role
    ) %>%
    filter(!is.na(myeloid_cluster)) %>%
    filter(myeloid_cluster %in% clusters_keep) %>%
    # Pair-level identity:
    mutate(pair_id = paste(ligand, receptor, sep = " | "))
}

subset_cluster <- function(df, cl) df %>% filter(myeloid_cluster == cl)

# Unique vs baseline set of pair IDs
unique_pairs_vs_baseline <- function(df_query, baseline_pair_set) {
  df_query %>% filter(!(pair_id %in% baseline_pair_set))
}

# Condense duplicates: keep max prob per (cluster, pair_id, Direction)
condense_unique <- function(df) {
  if (!("prob" %in% colnames(df))) {
    return(df %>% distinct(myeloid_cluster, pair_id, Direction, source, target, .keep_all = TRUE))
  }

  df %>%
    group_by(myeloid_cluster, pair_id, Direction) %>%
    arrange(desc(prob), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}

write_unique_interactions_workbook <- function(df_unique, out_xlsx) {
  wb <- createWorkbook()

  # Overview
  addWorksheet(wb, "All_batches")
  df_overview <- df_unique %>%
    transmute(
      MyeloidCluster = myeloid_cluster,
      Source = source,
      Target = target,
      Ligand = ligand,
      Receptor = receptor,
      PairID = pair_id,
      Direction, Role,
      Probability = prob
    ) %>%
    arrange(MyeloidCluster, Direction, PairID)
  writeData(wb, "All_batches", df_overview)

  # Per cluster sheets
  for (cl in clusters_keep) {
    df_cl <- subset_cluster(df_unique, cl)
    if (nrow(df_cl) == 0) next

    sheet_name <- paste0("batch_", cl)
    addWorksheet(wb, sheet_name)

    df_out <- df_cl %>%
      transmute(
        MyeloidCluster = myeloid_cluster,
        Source = source,
        Target = target,
        Ligand = ligand,
        Receptor = receptor,
        PairID = pair_id,
        Direction, Role,
        Probability = prob
      ) %>%
      arrange(Direction, PairID)

    writeData(wb, sheet = sheet_name, x = df_out)
  }

  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("Saved workbook: ", out_xlsx)
}

# Key used for comparing 1A-unique vs 1B-unique:
# include Direction so T->Myeloid vs Myeloid->T are treated as different.
make_key <- function(df) {
  df %>% mutate(key = paste(myeloid_cluster, pair_id, Direction, sep = " || "))
}

# ------------------------------------------------------------------------------
# 12D. Load and prepare data
# ------------------------------------------------------------------------------
df1A <- prep_df(read_all_sheet(comp1A_xlsx))
df1B <- prep_df(read_all_sheet(comp1B_xlsx))
df2 <- prep_df(read_all_sheet(comp2_xlsx))
df3 <- prep_df(read_all_sheet(comp3_xlsx))

# Baseline = union of Comp2 + Comp3 pair IDs (across kept clusters)
baseline_pairs_2_3 <- unique(c(df2$pair_id, df3$pair_id))

# ------------------------------------------------------------------------------
# 12E. Unique vs baseline
# ------------------------------------------------------------------------------
unique_1A <- df1A %>%
  unique_pairs_vs_baseline(baseline_pairs_2_3) %>%
  condense_unique()

unique_1B <- df1B %>%
  unique_pairs_vs_baseline(baseline_pairs_2_3) %>%
  condense_unique()

# Output files #1 and #2
out_unique_1A <- file.path(outdir_comp1A,
  "UNIQUE_interactions_Comp1A_vs_Comp2plusComp3_MonoMac_DC.xlsx"
)
out_unique_1B <- file.path(outdir_comp1B,
  "UNIQUE_interactions_Comp1B_vs_Comp2plusComp3_MonoMac_DC.xlsx"
)

write_unique_interactions_workbook(unique_1A, out_unique_1A)
write_unique_interactions_workbook(unique_1B, out_unique_1B)

# ------------------------------------------------------------------------------
# 12F. Third file: compare unique_1A vs unique_1B
# ------------------------------------------------------------------------------
u1A_flags <- make_key(unique_1A) %>%
  distinct(key) %>%
  mutate(In_1A_unique_vs_2plus3 = TRUE)

u1B_flags <- make_key(unique_1B) %>%
  distinct(key) %>%
  mutate(In_1B_unique_vs_2plus3 = TRUE)

membership <- full_join(u1A_flags, u1B_flags, by = "key") %>%
  mutate(
    In_1A_unique_vs_2plus3 = ifelse(is.na(In_1A_unique_vs_2plus3), FALSE, In_1A_unique_vs_2plus3),
    In_1B_unique_vs_2plus3 = ifelse(is.na(In_1B_unique_vs_2plus3), FALSE, In_1B_unique_vs_2plus3),
    UniquenessClass = case_when(
      In_1A_unique_vs_2plus3 & In_1B_unique_vs_2plus3 ~ "both_1A_and_1B_vs_2plus3",
      In_1A_unique_vs_2plus3 & !In_1B_unique_vs_2plus3 ~ "only_1A_vs_2plus3",
      !In_1A_unique_vs_2plus3 & In_1B_unique_vs_2plus3 ~ "only_1B_vs_2plus3",
      TRUE ~ NA_character_
    )
  )

combined_union <- bind_rows(
  make_key(unique_1A) %>% mutate(From = "Comp1A_unique_vs_2plus3"),
  make_key(unique_1B) %>% mutate(From = "Comp1B_unique_vs_2plus3")
) %>%
  # one representative row per key; keep highest prob if duplicates
  group_by(key) %>%
  arrange(desc(prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(membership %>% select(key, In_1A_unique_vs_2plus3, In_1B_unique_vs_2plus3, UniquenessClass),
            by = "key") %>%
  transmute(
    MyeloidCluster = myeloid_cluster,
    Source = source,
    Target = target,
    Ligand = ligand,
    Receptor = receptor,
    PairID = pair_id,
    Direction, Role,
    Probability = prob,
    In_1A_unique_vs_2plus3,
    In_1B_unique_vs_2plus3,
    UniquenessClass
  ) %>%
  arrange(MyeloidCluster, UniquenessClass, Direction, PairID)

out_unique_compare_1A_1B <- file.path(
  outdir_root,
  "UNIQUE_interactions_Comp1A_vs_Comp1B_classified_relative_to_Comp2plusComp3.xlsx"
)

wb3 <- createWorkbook()

# Overview
addWorksheet(wb3, "All_batches")
writeData(wb3, "All_batches", combined_union)

# Per cluster sheets
for (cl in clusters_keep) {
  df_cl <- combined_union %>% filter(MyeloidCluster == cl)
  if (nrow(df_cl) == 0) next
  sheet_name <- paste0("batch_", cl)
  addWorksheet(wb3, sheet_name)
  writeData(wb3, sheet_name, df_cl)
}

saveWorkbook(wb3, out_unique_compare_1A_1B, overwrite = TRUE)
message("Saved third comparison workbook: ", out_unique_compare_1A_1B)

