#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/visualization/biotin_violins_by_cluster.R
# Original file: 202060221 BIOTIN VIOLINS.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Generate per-cluster Biotin violin plots.
# Inputs: Canonical analysis-ready Seurat object with ADT, final_clusters, and condition metadata.
# Outputs: Per-cluster Biotin violin PDFs.
# Assay/layer input: ADT data for Biotin-TotalSeqC values.
# Dependencies: Seurat, ggplot2, dplyr.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Two plot variants are preserved: per-cluster violins and split-fill threshold violins.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!exists("obj")) {
  seurat_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
  seurat_rdata <- Sys.getenv("SEURAT_RDATA", unset = file.path(dirname(seurat_rds), "seurat_final_clean4.RData"))

  if (file.exists(seurat_rds)) {
    obj <- readRDS(seurat_rds)
  } else if (file.exists(seurat_rdata)) {
    load(seurat_rdata)
    if (!exists("obj") && exists("seurat_final_clean4")) obj <- seurat_final_clean4
  } else {
    stop("Could not find Seurat input. Set SEURAT_RDS or place seurat_final_clean4.rds at ", seurat_rds)
  }
}

stopifnot(exists("obj"))
stopifnot(inherits(obj, "Seurat"))

cluster_col <- "final_clusters"
cond_col <- "condition"
biotin_feat <- "Biotin-TotalSeqC"

treated_label <- "treated"
control_label <- "untreated"

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "20260221 BIOTIN VIOLINS_PER_CLUSTER")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Prepare metadata
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[cond_col]] <- as.character(obj[[cond_col]][,1])

clusters_plot <- sort(unique(obj[[cluster_col]][,1]))

# ------------------------------------------------------------------------------
# ADT CLR normalization (same as CellChat logic)
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)

adt <- GetAssayData(obj, assay = "ADT", slot = "data")
stopifnot(biotin_feat %in% rownames(adt))

# ------------------------------------------------------------------------------
# Build dataframe once
# ------------------------------------------------------------------------------
cells_use <- colnames(obj)[obj[[cond_col]][,1] %in% c(treated_label, control_label)]

df_all <- data.frame(
  cell = cells_use,
  cluster = obj[[cluster_col]][cells_use, 1],
  condition = obj[[cond_col]][cells_use, 1],
  biotin = as.numeric(adt[biotin_feat, cells_use, drop = TRUE]),
  stringsAsFactors = FALSE
)

df_all$condition_plot <- factor(
  ifelse(df_all$condition == treated_label, "Treated", "Control"),
  levels = c("Treated", "Control")
)

# ------------------------------------------------------------------------------
# Loop over clusters -> ONE PDF EACH
# ------------------------------------------------------------------------------
min_untreated_for_thr <- 5

for (cl in clusters_plot) {

  df <- df_all %>% filter(cluster == cl)

  # compute Q3 threshold from untreated cells
  untreated_vals <- df$biotin[df$condition == control_label]

  thr_q75 <- if (length(untreated_vals) >= min_untreated_for_thr) {
    as.numeric(quantile(untreated_vals, 0.75, na.rm = TRUE))
  } else {
    NA_real_
  }

  # plot
  p <- ggplot(df, aes(x = condition_plot, y = biotin, fill = condition_plot)) +
    geom_violin(trim = FALSE, width = 0.9, color = "grey30") +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      size = 0.6,
      alpha = 0.6,
      color = "grey20"
    ) +
    geom_hline(
      yintercept = thr_q75,
      linetype = "dotted",
      linewidth = 0.6,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "Treated" = "#7fb6e8",
        "Control" = "grey70"
      )
    ) +
    labs(
      title = paste0("Cluster ", cl, " - Biotin-TotalSeqC (ADT CLR)"),
      x = NULL,
      y = "Biotin-TotalSeqC (CLR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  outfile <- file.path(
    outdir,
    paste0("Violin_Biotin_CLR_cluster_", cl, "_treated_vs_control_Q3line.pdf")
  )

  pdf(outfile, width = 4.5, height = 5)
  print(p)
  dev.off()

  message("Saved cluster ", cl, ": ", outfile)
}

# restore default assay
DefaultAssay(obj) <- "RNA"

# color split by threshold

# ------------------------------------------------------------------------------
# 20260221 BIOTIN VIOLINS - treated violin split by threshold
# ------------------------------------------------------------------------------
# For EACH cluster (ONE PDF):
#   - LEFT:  Treated violin (single shape)
#       * below Q3 (control) = BLUE   (#859EE1)
#       * above Q3 (control) = YELLOW (#FFD84D)
#   - RIGHT: Control/Untreated violin = GREY (grey70)
#   - Dotted black line = Q3 threshold from untreated cells in that cluster
# ------------------------------------------------------------------------------
# Output folder:
#   output/20260221 BIOTIN VIOLINS_PER_CLUSTER_SPLITFILL
# ------------------------------------------------------------------------------

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
      stop("Loaded seurat_final_clean4.RData but neither 'obj' nor 'seurat_final_clean4' exists.")
    }
  }
} else {
  stop("Could not find Seurat input. Set SEURAT_RDS or place seurat_final_clean4.rds at ", seurat_rds)
}

stopifnot(inherits(obj, "Seurat"))
stopifnot(all(c("final_clusters", "condition") %in% colnames(obj@meta.data)))

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
cond_col <- "condition"
biotin_feat <- "Biotin-TotalSeqC"

treated_label <- "treated"
control_label <- "untreated"

min_untreated_for_thr <- 5

COL_BLUE <- "#859EE1"
COL_YELLOW <- "#FFD84D"
COL_GREY <- "grey70"
COL_EDGE <- "grey30"
COL_DOT <- "grey20"

outdir <- file.path(output_root, "20260221 BIOTIN VIOLINS_PER_CLUSTER_SPLITFILL")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1) PREP METADATA + ADT CLR
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][, 1])
obj[[cond_col]] <- as.character(obj[[cond_col]][, 1])

DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)

adt <- GetAssayData(obj, assay = "ADT", slot = "data")
stopifnot(biotin_feat %in% rownames(adt))

cells_use <- colnames(obj)[obj[[cond_col]][, 1] %in% c(treated_label, control_label)]

df_all <- data.frame(
  cell = cells_use,
  cluster = obj[[cluster_col]][cells_use, 1],
  condition = obj[[cond_col]][cells_use, 1],
  biotin = as.numeric(adt[biotin_feat, cells_use, drop = TRUE]),
  stringsAsFactors = FALSE
)

clusters_plot <- sort(unique(df_all$cluster))

# ------------------------------------------------------------------------------
# 2) HELPER: build a violin polygon and split at threshold
# ------------------------------------------------------------------------------
make_violin_split_polys <- function(values, x0, width = 0.45, thr) {
  # values: numeric vector
  # x0: x position
  # width: half-width (0.45 means total ~0.9)
  # thr: threshold where to split

  values <- values[is.finite(values)]
  if (length(values) < 2) return(list(below = NULL, above = NULL, full = NULL))

  d <- stats::density(values, n = 512, na.rm = TRUE)

  y <- d$x
  dens <- d$y

  # normalize to max density = 1 (like scale="width")
  dens_scaled <- dens / max(dens)
  halfw <- dens_scaled * width

  # Build full violin boundary
  x_right <- x0 + halfw
  x_left <- x0 - halfw

  full <- data.frame(
    x = c(x_right, rev(x_left)),
    y = c(y, rev(y)),
    part = "full"
  )

  # If thr is NA, return full only
  if (!is.finite(thr)) {
    return(list(below = NULL, above = NULL, full = full))
  }

  # helper: slice and close a polygon for a y-range
  build_part <- function(ymin, ymax, label) {
    idx <- which(y >= ymin & y <= ymax)
    if (length(idx) < 2) return(NULL)

    y_s <- y[idx]
    hr_s <- halfw[idx]

    # ensure boundary includes thr point if inside range
    add_thr <- function(vec_y, vec_hr) {
      if (thr > min(vec_y) && thr < max(vec_y)) {
        hr_thr <- approx(x = y, y = halfw, xout = thr)$y
        vec_y <- sort(c(vec_y, thr))
        # insert hr_thr aligned with thr
        # recompute hr values by interpolation for all vec_y
        vec_hr <- approx(x = y, y = halfw, xout = vec_y)$y
      }
      list(y = vec_y, hr = vec_hr)
    }

    tmp <- add_thr(y_s, hr_s)
    y_s <- tmp$y
    hr_s <- tmp$hr

    xr <- x0 + hr_s
    xl <- x0 - hr_s

    data.frame(
      x = c(xr, rev(xl)),
      y = c(y_s, rev(y_s)),
      part = label
    )
  }

  below <- build_part(min(y), min(thr, max(y)), "below")
  above <- build_part(max(min(thr, max(y)), min(y)), max(y), "above")

  list(below = below, above = above, full = full)
}

# ------------------------------------------------------------------------------
# 3) LOOP: ONE PDF PER CLUSTER
# ------------------------------------------------------------------------------
for (cl in clusters_plot) {

  df <- df_all %>% filter(cluster == cl)

  df_treat <- df %>% filter(condition == treated_label)
  df_ctrl <- df %>% filter(condition == control_label)

  # Q3 threshold from untreated cells
  thr_q75 <- if (nrow(df_ctrl) >= min_untreated_for_thr) {
    as.numeric(quantile(df_ctrl$biotin, 0.75, na.rm = TRUE))
  } else {
    NA_real_
  }

  # Build polygons
  # x positions: Treated = 1, Control = 2
  polys_t <- make_violin_split_polys(df_treat$biotin, x0 = 1, width = 0.45, thr = thr_q75)
  polys_c <- make_violin_split_polys(df_ctrl$biotin,  x0 = 2, width = 0.45, thr = NA_real_) # no split

  # Dots: color treated by above/below threshold; control grey
  df_pts <- df %>%
    mutate(
      x = ifelse(condition == treated_label, 1, 2),
      dot_group = case_when(
        condition == control_label ~ "Control",
        is.na(thr_q75) ~ "Treated <= Q3",
        biotin > thr_q75 ~ "Treated > Q3",
        TRUE ~ "Treated <= Q3"
      )
    )

  # Base plot
  p <- ggplot() +
    # Control violin (grey)
    {if (!is.null(polys_c$full)) geom_polygon(
      data = polys_c$full,
      aes(x = x, y = y),
      fill = COL_GREY,
      color = COL_EDGE,
      linewidth = 0.3
    )} +

    # Treated violin: below (blue) + above (yellow)
    {if (!is.null(polys_t$below)) geom_polygon(
      data = polys_t$below,
      aes(x = x, y = y),
      fill = COL_BLUE,
      color = NA
    )} +
    {if (!is.null(polys_t$above)) geom_polygon(
      data = polys_t$above,
      aes(x = x, y = y),
      fill = COL_YELLOW,
      color = NA
    )} +
    # Outline treated violin once (so it looks like ONE violin)
    {if (!is.null(polys_t$full)) geom_polygon(
      data = polys_t$full,
      aes(x = x, y = y),
      fill = NA,
      color = COL_EDGE,
      linewidth = 0.3
    )} +

    # Points
    geom_point(
      data = df_pts,
      aes(
        x = x,
        y = biotin,
        color = dot_group
      ),
      position = position_jitter(width = 0.10, height = 0),
      size = 0.6,
      alpha = 0.6
    ) +

    # Threshold line
    geom_hline(
      yintercept = thr_q75,
      linetype = "dotted",
      linewidth = 0.6,
      color = "black"
    ) +

    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Treated", "Control"),
      limits = c(0.4, 2.6)
    ) +

    scale_color_manual(
      values = c(
        "Treated <= Q3" = COL_DOT,
        "Treated > Q3" = COL_DOT,
        "Control" = COL_DOT
      )
    ) +

    labs(
      title = paste0("Cluster ", cl, " - Biotin-TotalSeqC (ADT CLR)"),
      subtitle = ifelse(
        is.na(thr_q75),
        "Q3 threshold: NA (too few untreated cells)",
        paste0("Q3 threshold (untreated): ", signif(thr_q75, 4))
      ),
      x = NULL,
      y = "Biotin-TotalSeqC (CLR)"
    ) +

    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )

  outfile <- file.path(outdir, paste0("Violin_Biotin_CLR_cluster_", cl, "_splitFillTreated.pdf"))

  pdf(outfile, width = 4.5, height = 5)
  print(p)
  dev.off()

  message("Saved cluster ", cl, ": ", outfile)
}

DefaultAssay(obj) <- "RNA"
message("Done. Output folder: ", outdir)

