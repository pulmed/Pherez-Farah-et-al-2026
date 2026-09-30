#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/visualization/plot_pseudobulk_volcano_panels.R
# Original file: 20260218 VOLCANOS.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Replot selected pseudobulk edgeR comparisons as condensed volcano panels.
# Inputs: edgeR XLSX outputs from scripts/pseudobulk/run_cluster_pseudobulk_dge.R.
# Outputs: Individual and combined volcano-panel PDFs/PNGs.
# Assay/layer input: Uses precomputed RNA pseudobulk DGE result tables.
# Dependencies: readxl, ggplot2, ggrepel, patchwork.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - This script preserves the final display settings from the selected manuscript volcano panels.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

dge_dir <- Sys.getenv("DGE_DIR", unset = file.path("output", "16_DGE_clusters_RNA_pseudobulk"))
output_dir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "15_pseudobulk_volcano_panels")

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !dir.exists(dge_dir)) {
  message("SMOKE_TEST=1: DGE directory not found; skipping volcano-panel plotting.")
  quit(save = "no", status = 0)
}

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: DGE directory found: ", dge_dir)
  message("SMOKE_TEST=1: skipping volcano-panel plotting.")
  quit(save = "no", status = 0)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fc_threshold <- 1
fdr_threshold <- 0.01
y_threshold <- -log10(fdr_threshold)
x_limit <- 4

colors <- c(
  ns = "#2b6cb0",
  sig = "#c53030",
  highlight = "#2f6f55"
)

comparisons <- list(
  list(
    name = "Cluster2_vs_Cluster4",
    title = "Cluster2 vs Cluster4",
    file = file.path(dge_dir, "Cluster2_vs_Cluster4_edgeR.xlsx"),
    y_cap = 20,
    genes = c(
      "Ccl8", "RetnIa", "Folr2", "C1qa", "Ccl7", "C1qc", "C1qb",
      "Plet1", "Ly6a2", "Il1r2", "Sell", "Trem1", "Plac8", "Ly6c2",
      "Ccr2", "Prtn3", "Fabp4", "Arg1", "Spp1", "Fabp5", "Hpgds",
      "Mmp12", "Ccl5", "Cd36", "Apoe", "Trem2", "Mertk", "Ccl12",
      "Ccl24", "Pf4", "Flrt2", "Pdpn"
    )
  ),
  list(
    name = "Cluster0_vs_Cluster5",
    title = "Cluster0 vs Cluster5",
    file = file.path(dge_dir, "Cluster0_vs_Cluster5_edgeR.xlsx"),
    y_cap = 12,
    genes = c(
      "Il1b", "Ccl5", "Ly6a2", "Ly6c2", "Ly6i", "Ly6a", "Cxcl9",
      "Isoc1", "Plac8", "IL18", "Nlrp1b", "Ly6e", "Tap2Il1b",
      "Il1a", "Tnf", "Nos2", "Il6", "Ccl7", "Fabp4", "Cd200",
      "Spp1", "Kdr", "Fabp5", "Ccl9", "Dhfr", "Ccl6", "Ccl12",
      "Ccl24", "Pf4", "Flrt2", "Cd36", "Pdpn", "Mmp12", "Ccl2",
      "Ccl8", "Ccl4"
    )
  )
)

detect_column <- function(data, candidates, label) {
  hit <- candidates[candidates %in% colnames(data)][1]
  if (is.na(hit)) stop("Could not detect ", label, " column in: ", paste(colnames(data), collapse = ", "))
  hit
}

make_volcano <- function(config) {
  if (!file.exists(config$file)) stop("Missing DGE workbook: ", config$file)

  df <- readxl::read_excel(config$file)
  if (!"gene" %in% colnames(df)) df$gene <- as.character(df[[1]])

  fc_col <- detect_column(df, c("log2FC", "avg_log2FC", "logFC", "log2FoldChange"), "fold-change")
  fdr_col <- detect_column(df, c("FDR", "padj", "p_val_adj", "adj.P.Val"), "FDR")

  df$log2FC <- as.numeric(df[[fc_col]])
  df$FDR <- as.numeric(df[[fdr_col]])
  df$neglog10FDR <- -log10(pmax(df$FDR, .Machine$double.xmin))
  df$neglog10FDR[df$neglog10FDR > config$y_cap] <- config$y_cap

  highlight_genes <- unique(config$genes)
  df$group <- "ns"
  df$group[df$FDR <= fdr_threshold & abs(df$log2FC) >= fc_threshold] <- "sig"
  df$group[df$gene %in% highlight_genes] <- "highlight"
  label_df <- df[df$gene %in% highlight_genes, , drop = FALSE]

  ggplot(df, aes(log2FC, neglog10FDR)) +
    geom_point(aes(color = group), size = 2.2, alpha = 0.9, na.rm = TRUE) +
    geom_point(
      data = df[df$group == "highlight", , drop = FALSE],
      color = colors["highlight"],
      size = 3.6,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dotted", linewidth = 0.4) +
    geom_hline(yintercept = y_threshold, linetype = "dotted", linewidth = 0.4) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 5,
      color = "black",
      box.padding = 0.5,
      point.padding = 0.35,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 1
    ) +
    scale_color_manual(values = colors) +
    coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, config$y_cap), expand = FALSE) +
    labs(title = config$title, x = "log2 Fold Change", y = expression(-log[10]("FDR"))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.6)
    )
}

plots <- lapply(comparisons, make_volcano)
names(plots) <- vapply(comparisons, function(x) x$name, character(1))

for (name in names(plots)) {
  ggsave(file.path(output_dir, paste0("Volcano_", name, ".pdf")), plots[[name]], width = 5.2, height = 5.2)
  ggsave(file.path(output_dir, paste0("Volcano_", name, ".png")), plots[[name]], width = 5.2, height = 5.2, dpi = 400)
}

combined <- plots[[1]] + plots[[2]] + plot_layout(ncol = 2)
ggsave(file.path(output_dir, "Figure_Volcanos_2panel.pdf"), combined, width = 10.4, height = 5.2, useDingbats = FALSE)
ggsave(file.path(output_dir, "Figure_Volcanos_2panel.png"), combined, width = 10.4, height = 5.2, dpi = 400)

message("Done. Volcano panels written to: ", output_dir)
