#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_generate_general_metrics_overview.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Generate general overview metrics, plots, marker summaries, heatmaps, and selected signature visualizations.
# Inputs: Normalized Seurat object from 02_demultiplex_normalize_and_visualize_seurat.R.
# Outputs: PNG plots and tabular summaries.
# Dependencies: Seurat, dplyr, tidyr, ggplot2, scales, openxlsx, pheatmap.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Edit the configuration section before running.
# - RNA log-normalization is expected to have been performed in script 02.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(openxlsx)
  library(pheatmap)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION  ← EDIT THIS SECTION
# ------------------------------------------------------------------------------

input_file <- "path/to/combined_seurat_normalized.rds"
output_dir <- "path/to/general_metrics_overview"

cluster_column <- "final_clusters"
condition_column <- "condition"

rna_assay <- "RNA"
sct_assay <- "SCT"
umap_reduction <- "umap"

condition_colors <- c(
  treated = "#7C9FE6",
  untreated = "gray80",
  control = "gray80"
)

cluster_order <- c(
  "0", "2", "3", "4", "5", "14", "17",
  "13", "1", "7", "8", "9",
  "6", "10", "11", "15",
  "12", "16", "18", "19"
)

features_to_plot <- c("Cd36", "H2-Ab1", "Adgre1", "Nos2", "Arg1", "Cxcl16", "Ccr5")
features_split_by_condition <- c("Cxcl16", "Ccr5")

marker_gene_for_stats <- "Ly6c2"
marker_stats_assays <- c("RNA", "SCT")

run_top_markers <- TRUE
top_marker_cluster <- "7"
top_marker_n <- 50
top_marker_min_pct <- 0.25
top_marker_logfc_threshold <- 0.25

run_marker_heatmap <- TRUE
heatmap_assay <- "RNA"
heatmap_clusters <- c("0", "2", "3", "4", "5", "13", "14", "17")
heatmap_genes <- c(
  "Nos2", "Slc7a2", "Arg1", "Mgst1", "Fpr2", "Il1a", "Slpi", "Inhba",
  "Grina", "Ly6i", "Cxcl2", "Prdx1", "Cfb", "Pla2g7", "Blvrb", "Ptgs2",
  "Slc7a11", "Cacna1d", "Trem1", "Arg2", "F10", "Srxn1", "Anpep", "Bst1",
  "Tarm1", "Fbxl5", "Prdx5", "Tnf", "Carmil1", "Clec4d", "Acod1",
  "Clec4e", "Tgfbi", "Smpdl3b", "Sod2", "Sgk1", "Hcar2", "Upp1",
  "Bnip3", "Msr1", "Met", "Sgms2", "Pla2g4a", "Prdx6", "Slc7a8",
  "Fth1", "Gatm", "Glrx", "Ccrl2", "Hk3", "Nlrp3", "Cd14", "Cpd",
  "Acsl1", "Pstpip2", "Gclm", "Fmnl2", "Cyp4f18", "Ly6c2"
)

run_three_signature_overlay <- TRUE
cdc1_genes <- c("Xcr1", "Clec9a", "Itgae", "Cadm1")
cdc2_genes <- c("Cd209a", "Clec10a", "Mgl2")
mature_dc_genes <- c("Ccr7", "Fscn1", "Il12b", "Cacnb3")
signature_z_cutoff <- 1.0
signature_quantile <- 0.90

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

save_png <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

has_assay <- function(object, assay_name) {
  assay_name %in% Assays(object)
}

require_metadata <- function(object, columns) {
  missing_columns <- setdiff(columns, colnames(object@meta.data))
  if (length(missing_columns) > 0) {
    stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
  }
}

set_cluster_idents <- function(object, cluster_column) {
  object[[cluster_column]] <- as.character(object[[cluster_column]][, 1])
  Idents(object) <- object[[cluster_column]][, 1]
  object
}

clean_condition_labels <- function(x) {
  x <- tolower(trimws(as.character(x)))

  dplyr::recode(
    x,
    "ctrl" = "control",
    "control" = "control",
    "untreated" = "control",
    "vehicle" = "control",
    "mock" = "control",
    "treated" = "treated",
    "treat" = "treated",
    "drug" = "treated",
    .default = x
  )
}

write_xlsx <- function(tables, filename) {
  wb <- createWorkbook()

  for (sheet_name in names(tables)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, tables[[sheet_name]])
  }

  saveWorkbook(wb, filename, overwrite = TRUE)
}

get_umap_dataframe <- function(object, reduction, cluster_column, condition_column) {
  umap_df <- as.data.frame(Embeddings(object, reduction))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

  umap_df$cell <- rownames(umap_df)
  umap_df$cluster <- as.character(object@meta.data[umap_df$cell, cluster_column])
  umap_df$condition <- as.character(object@meta.data[umap_df$cell, condition_column])

  umap_df
}

# ------------------------------------------------------------------------------
# 3. LOAD DATA AND VALIDATE INPUTS
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

output_dir <- make_dir(output_dir)

seu <- readRDS(input_file)

require_metadata(seu, c(cluster_column, condition_column))

if (!umap_reduction %in% Reductions(seu)) {
  stop("UMAP reduction not found: ", umap_reduction)
}

seu <- set_cluster_idents(seu, cluster_column)

message("Loaded Seurat object with ", ncol(seu), " cells.")
message("Assays: ", paste(Assays(seu), collapse = ", "))
message("Clusters: ", paste(sort(unique(Idents(seu))), collapse = ", "))

# ------------------------------------------------------------------------------
# 4. UMAP OVERVIEW PLOTS
# ------------------------------------------------------------------------------

umap_dir <- make_dir(file.path(output_dir, "umap_overview"))

p_cluster <- DimPlot(
  seu,
  reduction = umap_reduction,
  group.by = cluster_column,
  label = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  ggtitle("UMAP by Cluster")

p_condition <- DimPlot(
  seu,
  reduction = umap_reduction,
  group.by = condition_column,
  pt.size = 0.5
) +
  scale_color_manual(values = condition_colors, na.value = "gray80") +
  ggtitle("UMAP by Condition") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_split <- DimPlot(
  seu,
  reduction = umap_reduction,
  split.by = condition_column,
  group.by = cluster_column,
  pt.size = 0.5
) +
  ggtitle("UMAP Split by Condition")

save_png(p_cluster, file.path(umap_dir, "UMAP_by_cluster.png"))
save_png(p_condition, file.path(umap_dir, "UMAP_by_condition.png"))
save_png(p_split, file.path(umap_dir, "UMAP_split_by_condition.png"), width = 10)

# ------------------------------------------------------------------------------
# 5. CLUSTER HIGHLIGHT PLOTS
# ------------------------------------------------------------------------------

highlight_dir <- make_dir(file.path(output_dir, "cluster_condition_highlights"))

umap_df <- get_umap_dataframe(
  object = seu,
  reduction = umap_reduction,
  cluster_column = cluster_column,
  condition_column = condition_column
)

clusters <- sort(unique(umap_df$cluster))

for (cluster_id in clusters) {
  plot_df <- umap_df %>%
    mutate(
      highlight_condition = ifelse(cluster == cluster_id, condition, "Other"),
      highlight_condition = factor(
        highlight_condition,
        levels = c("Other", "treated", "untreated", "control")
      )
    )

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = highlight_condition)) +
    geom_point(size = 0.4, alpha = 0.85) +
    scale_color_manual(
      values = c(
        Other = "gray85",
        treated = "#E76F51",
        untreated = "#2A9D8F",
        control = "gray60"
      ),
      drop = FALSE
    ) +
    coord_equal() +
    theme_void(base_size = 14) +
    ggtitle(paste("Cluster", cluster_id, "- Condition Overlay"))

  save_png(
    p,
    file.path(highlight_dir, paste0("Cluster_", cluster_id, "_condition_highlight.png"))
  )
}

# ------------------------------------------------------------------------------
# 6. CLUSTER COMPOSITION AND ENRICHMENT
# ------------------------------------------------------------------------------

composition_dir <- make_dir(file.path(output_dir, "cluster_condition_composition"))

metadata <- as.data.frame(seu@meta.data) %>%
  mutate(
    cluster = as.character(.data[[cluster_column]]),
    condition_raw = as.character(.data[[condition_column]]),
    condition_clean = clean_condition_labels(condition_raw)
  )

cluster_counts <- metadata %>%
  count(cluster, condition_clean, name = "n")

cluster_props_within_cluster <- cluster_counts %>%
  group_by(cluster) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

cluster_levels <- intersect(cluster_order, unique(cluster_props_within_cluster$cluster))

if (length(cluster_levels) > 0) {
  cluster_props_within_cluster$cluster <- factor(
    cluster_props_within_cluster$cluster,
    levels = cluster_levels
  )
}

p_composition <- ggplot(
  cluster_props_within_cluster,
  aes(x = cluster, y = fraction, fill = condition_clean)
) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(
    values = c(control = "#EDEDED", treated = "#AAB8FF"),
    drop = FALSE,
    name = "Condition"
  ) +
  labs(
    title = "Condition Composition Within Each Cluster",
    x = "Cluster",
    y = "Fraction of Cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_png(
  p_composition,
  file.path(composition_dir, "cluster_condition_composition.png"),
  width = 8
)

cluster_props_within_condition <- cluster_counts %>%
  group_by(condition_clean) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

props_wide <- cluster_props_within_condition %>%
  select(cluster, condition_clean, fraction) %>%
  pivot_wider(
    names_from = condition_clean,
    values_from = fraction,
    values_fill = 0
  )

if (all(c("treated", "control") %in% colnames(props_wide))) {
  enrichment <- props_wide %>%
    mutate(
      log2_enrichment = log2((treated + 1e-6) / (control + 1e-6)),
      enriched_group = ifelse(log2_enrichment > 0, "treated", "control")
    ) %>%
    arrange(desc(log2_enrichment))

  enrichment$cluster <- factor(enrichment$cluster, levels = enrichment$cluster)

  p_enrichment <- ggplot(
    enrichment,
    aes(x = cluster, y = log2_enrichment, fill = enriched_group)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(
      values = c(treated = "#7C9FE6", control = "gray80"),
      name = NULL
    ) +
    labs(
      title = "Cluster Enrichment by Condition",
      x = "Cluster",
      y = "log2(treated / control)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )

  save_png(
    p_enrichment,
    file.path(composition_dir, "cluster_log2_enrichment.png")
  )

  write.csv(
    enrichment,
    file.path(composition_dir, "cluster_log2_enrichment.csv"),
    row.names = FALSE
  )
}

write.csv(
  cluster_counts,
  file.path(composition_dir, "cluster_condition_counts.csv"),
  row.names = FALSE
)

write.csv(
  cluster_props_within_condition,
  file.path(composition_dir, "cluster_condition_fractions_within_condition.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. FEATURE AND VIOLIN PLOTS
# ------------------------------------------------------------------------------

feature_dir <- make_dir(file.path(output_dir, "feature_visualization"))

if (has_assay(seu, rna_assay)) {
  DefaultAssay(seu) <- rna_assay

  available_features <- intersect(features_to_plot, rownames(seu[[rna_assay]]))

  for (feature in available_features) {
    p_feature <- FeaturePlot(
      seu,
      features = feature,
      reduction = umap_reduction,
      max.cutoff = "q95"
    ) +
      ggtitle(feature)

    p_vln_cluster <- VlnPlot(
      seu,
      features = feature,
      group.by = cluster_column,
      pt.size = 0
    ) +
      ggtitle(paste(feature, "by Cluster"))

    save_png(
      p_feature,
      file.path(feature_dir, paste0(feature, "_featureplot.png"))
    )

    save_png(
      p_vln_cluster,
      file.path(feature_dir, paste0(feature, "_violin_by_cluster.png")),
      width = 9
    )

    if (feature %in% features_split_by_condition) {
      p_split_feature <- FeaturePlot(
        seu,
        features = feature,
        reduction = umap_reduction,
        split.by = condition_column,
        max.cutoff = "q95"
      )

      p_vln_condition <- VlnPlot(
        seu,
        features = feature,
        group.by = condition_column,
        pt.size = 0
      ) +
        ggtitle(paste(feature, "by Condition"))

      save_png(
        p_split_feature,
        file.path(feature_dir, paste0(feature, "_featureplot_split_by_condition.png")),
        width = 10
      )

      save_png(
        p_vln_condition,
        file.path(feature_dir, paste0(feature, "_violin_by_condition.png"))
      )
    }
  }
}

# ------------------------------------------------------------------------------
# 8. MARKER STATISTICS
# ------------------------------------------------------------------------------

stats_dir <- make_dir(file.path(output_dir, paste0(marker_gene_for_stats, "_cluster_stats")))

run_marker_stats <- function(seu, assay, gene, cluster_column, out_dir) {
  if (!has_assay(seu, assay)) {
    message("Skipping marker stats; assay not found: ", assay)
    return(NULL)
  }

  if (!gene %in% rownames(seu[[assay]])) {
    message("Skipping marker stats; gene not found in ", assay, ": ", gene)
    return(NULL)
  }

  DefaultAssay(seu) <- assay

  df <- FetchData(seu, vars = c(gene, cluster_column)) %>%
    rename(expr = all_of(gene), cluster = all_of(cluster_column)) %>%
    mutate(cluster = factor(cluster))

  summary_table <- df %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean = mean(expr, na.rm = TRUE),
      median = median(expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median))

  kw <- kruskal.test(expr ~ cluster, data = df)

  kw_table <- data.frame(
    assay = assay,
    gene = gene,
    statistic = unname(kw$statistic),
    parameter = unname(kw$parameter),
    p_value = kw$p.value,
    method = kw$method
  )

  write.csv(
    summary_table,
    file.path(out_dir, paste0(gene, "_summary_by_cluster_", assay, ".csv")),
    row.names = FALSE
  )

  write.csv(
    kw_table,
    file.path(out_dir, paste0(gene, "_kruskal_wallis_", assay, ".csv")),
    row.names = FALSE
  )

  write_xlsx(
    list(
      Summary_by_cluster = summary_table,
      Kruskal_Wallis = kw_table
    ),
    file.path(out_dir, paste0(gene, "_stats_", assay, ".xlsx"))
  )

  invisible(list(summary = summary_table, kruskal = kw_table))
}

for (assay in marker_stats_assays) {
  run_marker_stats(
    seu = seu,
    assay = assay,
    gene = marker_gene_for_stats,
    cluster_column = cluster_column,
    out_dir = stats_dir
  )
}

# ------------------------------------------------------------------------------
# 9. TOP MARKERS FOR SELECTED CLUSTER
# ------------------------------------------------------------------------------

if (run_top_markers && has_assay(seu, sct_assay)) {
  marker_dir <- make_dir(file.path(output_dir, "top_cluster_markers"))

  DefaultAssay(seu) <- sct_assay
  seu <- set_cluster_idents(seu, cluster_column)

  if (top_marker_cluster %in% levels(Idents(seu))) {
    markers <- FindMarkers(
      seu,
      ident.1 = top_marker_cluster,
      min.pct = top_marker_min_pct,
      logfc.threshold = top_marker_logfc_threshold,
      only.pos = TRUE
    ) %>%
      tibble::rownames_to_column("gene") %>%
      filter(
        !grepl("^mt-", gene, ignore.case = TRUE),
        !grepl("^Rpl|^Rps", gene, ignore.case = TRUE)
      ) %>%
      arrange(desc(avg_log2FC))

    top_markers <- head(markers, top_marker_n)

    write.csv(
      top_markers,
      file.path(marker_dir, paste0("top", top_marker_n, "_cluster_", top_marker_cluster, "_markers.csv")),
      row.names = FALSE
    )

    write_xlsx(
      list(Top_markers = top_markers),
      file.path(marker_dir, paste0("top", top_marker_n, "_cluster_", top_marker_cluster, "_markers.xlsx"))
    )
  } else {
    message("Top-marker cluster not found: ", top_marker_cluster)
  }
}

# ------------------------------------------------------------------------------
# 10. SELECTED GENE HEATMAP
# ------------------------------------------------------------------------------

if (run_marker_heatmap && has_assay(seu, heatmap_assay)) {
  heatmap_dir <- make_dir(file.path(output_dir, "selected_gene_heatmap"))

  DefaultAssay(seu) <- heatmap_assay
  seu <- set_cluster_idents(seu, cluster_column)

  keep_clusters <- intersect(heatmap_clusters, levels(Idents(seu)))
  genes_present <- intersect(heatmap_genes, rownames(seu[[heatmap_assay]]))

  if (length(keep_clusters) == 0) {
    message("Skipping heatmap; none of the requested clusters are present.")
  } else if (length(genes_present) == 0) {
    message("Skipping heatmap; none of the requested genes are present.")
  } else {
    seu_sub <- subset(seu, idents = keep_clusters)

    avg_list <- AverageExpression(
      seu_sub,
      assays = heatmap_assay,
      features = genes_present,
      slot = "data",
      group.by = cluster_column
    )

    avg <- avg_list[[heatmap_assay]]

    col_ids <- sub("^.*?(-?\\d+)$", "\\1", colnames(avg))
    col_map <- setNames(colnames(avg), col_ids)
    present_cols <- unname(col_map[intersect(heatmap_clusters, names(col_map))])
    present_cols <- present_cols[!is.na(present_cols)]

    avg <- avg[, present_cols, drop = FALSE]

    z_scores <- t(scale(t(as.matrix(avg))))
    z_scores[is.na(z_scores)] <- 0

    palette <- colorRampPalette(c("#2b5cbe", "white", "#b40426"))(100)

    png(
      file.path(heatmap_dir, paste0("selected_genes_heatmap_", heatmap_assay, ".png")),
      width = 1400,
      height = 2000,
      res = 200
    )
    pheatmap(
      z_scores,
      color = palette,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      fontsize_row = 7,
      fontsize_col = 10,
      main = paste0("Selected genes across clusters – ", heatmap_assay),
      angle_col = 45
    )
    dev.off()

    write_xlsx(
      list(
        Average_expression = data.frame(Gene = rownames(avg), avg, row.names = NULL),
        Row_z_scores = data.frame(Gene = rownames(z_scores), z_scores, row.names = NULL)
      ),
      file.path(heatmap_dir, paste0("selected_genes_heatmap_", heatmap_assay, ".xlsx"))
    )
  }
}

# ------------------------------------------------------------------------------
# 11. SIGNATURE OVERLAY
# ------------------------------------------------------------------------------

if (run_three_signature_overlay && has_assay(seu, rna_assay)) {
  signature_dir <- make_dir(file.path(output_dir, "three_signature_overlay"))

  DefaultAssay(seu) <- rna_assay

  rna_data <- GetAssayData(seu, assay = rna_assay, slot = "data")

  score_signature <- function(genes) {
    genes_present <- intersect(genes, rownames(rna_data))

    if (length(genes_present) == 0) {
      return(rep(FALSE, ncol(seu)))
    }

    matrix <- t(as.matrix(rna_data[genes_present, , drop = FALSE]))
    z <- scale(matrix)
    z[is.na(z)] <- 0

    per_gene_positive <- z >= signature_z_cutoff
    hard_positive <- rowSums(per_gene_positive) >= min(2, length(genes_present))

    soft_score <- rowMeans(pmax(z, 0))
    soft_positive <- soft_score >= quantile(soft_score, signature_quantile, na.rm = TRUE)

    hard_positive & soft_positive
  }

  high_cdc1 <- score_signature(cdc1_genes)
  high_cdc2 <- score_signature(cdc2_genes)
  high_mature <- score_signature(mature_dc_genes)

  names(high_cdc1) <- colnames(seu)
  names(high_cdc2) <- colnames(seu)
  names(high_mature) <- colnames(seu)

  plot_df <- get_umap_dataframe(
    object = seu,
    reduction = umap_reduction,
    cluster_column = cluster_column,
    condition_column = condition_column
  )

  plot_df$cdc1 <- high_cdc1[plot_df$cell]
  plot_df$cdc2 <- high_cdc2[plot_df$cell]
  plot_df$mature <- high_mature[plot_df$cell]

  plot_df$state <- case_when(
    plot_df$cdc1 & !plot_df$cdc2 & !plot_df$mature ~ "cDC1 only",
    !plot_df$cdc1 & plot_df$cdc2 & !plot_df$mature ~ "cDC2 only",
    !plot_df$cdc1 & !plot_df$cdc2 & plot_df$mature ~ "Mature only",
    plot_df$cdc1 & plot_df$cdc2 & !plot_df$mature ~ "cDC1+cDC2",
    plot_df$cdc1 & !plot_df$cdc2 & plot_df$mature ~ "cDC1+Mature",
    !plot_df$cdc1 & plot_df$cdc2 & plot_df$mature ~ "cDC2+Mature",
    plot_df$cdc1 & plot_df$cdc2 & plot_df$mature ~ "All three",
    TRUE ~ "Other"
  )

  plot_df$state <- factor(
    plot_df$state,
    levels = c(
      "Other",
      "cDC1 only",
      "cDC2 only",
      "Mature only",
      "cDC1+cDC2",
      "cDC1+Mature",
      "cDC2+Mature",
      "All three"
    )
  )

  signature_colors <- c(
    "Other" = "#E6E9EF",
    "cDC1 only" = "#7AA874",
    "cDC2 only" = "#F4C095",
    "Mature only" = "#5C6AC4",
    "cDC1+cDC2" = "#B3C18C",
    "cDC1+Mature" = "#6D8FB6",
    "cDC2+Mature" = "#9FA2CF",
    "All three" = "#4B3F72"
  )

  p_signature <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = state)) +
    geom_point(size = 0.4, alpha = 0.85) +
    scale_color_manual(values = signature_colors, drop = FALSE) +
    coord_equal() +
    theme_void(base_size = 14) +
    ggtitle("UMAP RNA Signature Overlay: cDC1 / cDC2 / Mature DC")

  save_png(
    p_signature,
    file.path(signature_dir, "three_signature_umap_overlay.png"),
    width = 8
  )

  write.csv(
    as.data.frame(table(plot_df$state)),
    file.path(signature_dir, "three_signature_state_counts.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 12. FINISH
# ------------------------------------------------------------------------------

message("Done. Outputs saved to: ", output_dir)