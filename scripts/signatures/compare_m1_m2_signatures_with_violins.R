#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/compare_m1_m2_signatures_with_violins.R
# Original file: 20260824 M1 M2 FINAL.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Compare final M1-like and M2-like macrophage signatures in clusters 0 and 5.
# Inputs: Canonical analysis-ready Seurat object with RNA data and cluster/sample metadata.
# Outputs: Violin plots, paired sample-level plots, CSV summaries, and Excel workbooks.
# Assay/layer input: RNA data.
# Dependencies: Seurat, SeuratObject, dplyr, tidyr, ggplot2, ggridges, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Preserves legacy M1_M2_cluster0_vs5 output names.
# - Uses final 2026-08-24 M1/M2 gene sets.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(openxlsx)
})

source("scripts/utils/seurat_io.R")
source("scripts/utils/signature_helpers.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "14_M1_VS_M2", "M1_M2_final"))

assay_use <- Sys.getenv("ASSAY", unset = "RNA")
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "final_clusters")
sample_col_candidates <- parse_csv_env("SAMPLE_COL_CANDIDATES", c("HTO_classification", "hash.ID"))
clusters_use <- parse_csv_env("CLUSTERS_USE", c("0", "5"))
selected_umap_clusters <- parse_csv_env("SELECTED_UMAP_CLUSTERS", c("17", "0", "5", "2", "14", "3", "4"))

m1_genes <- c("Nos2", "Il1a", "Il1b", "Tnf", "Il6", "Cd80", "Cd86", "Cxcl9", "Cxcl10", "Stat1")
m2_genes <- c("Arg1", "Mrc1", "Retnla", "Chil3", "Il10", "Tgfb1", "Cd163", "Ccl17", "Ccl22", "Stat6")

score_seed <- as.integer(Sys.getenv("MODULE_SCORE_SEED", unset = "1234"))
save_scored_rds <- identical(tolower(Sys.getenv("SAVE_SCORED_RDS", unset = "false")), "true")

# ------------------------------------------------------------------------------
# LOAD AND VALIDATE
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_assays(obj, assay_use)
require_metadata(obj, cluster_col)
sample_col <- first_existing_metadata(obj, sample_col_candidates)

DefaultAssay(obj) <- assay_use
obj <- coerce_metadata_character(obj, c(cluster_col, sample_col))
obj <- ensure_rna_data_layer(obj, assay = assay_use)

features <- rownames(obj[[assay_use]])
m1_present <- present_signature_genes(m1_genes, features, "M1")
m2_present <- present_signature_genes(m2_genes, features, "M2")

gene_check <- bind_rows(
  gene_presence_table(m1_genes, features, "M1"),
  gene_presence_table(m2_genes, features, "M2")
)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: M1 genes present: ", length(m1_present), "/", length(m1_genes))
  message("SMOKE_TEST=1: M2 genes present: ", length(m2_present), "/", length(m2_genes))
  message("SMOKE_TEST=1: sample column selected: ", sample_col)
  quit(save = "no", status = 0)
}

write.csv(
  gene_check,
  file.path(outdir, "M1_M2_signature_gene_presence.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# MODULE SCORES
# ------------------------------------------------------------------------------
m1_score <- add_module_score(obj, m1_present, "M1_Final", assay = assay_use, seed = score_seed)
obj <- m1_score$object
m2_score <- add_module_score(obj, m2_present, "M2_Final", assay = assay_use, seed = score_seed)
obj <- m2_score$object

# ------------------------------------------------------------------------------
# CELL-LEVEL DATA
# ------------------------------------------------------------------------------
Idents(obj) <- obj[[cluster_col]][, 1]
obj_05 <- subset(obj, idents = clusters_use)

df_cell_wide <- FetchData(
  obj_05,
  vars = c(m1_score$column, m2_score$column, cluster_col, sample_col)
) %>%
  rename(
    M1_score = all_of(m1_score$column),
    M2_score = all_of(m2_score$column),
    cluster = all_of(cluster_col),
    sample = all_of(sample_col)
  ) %>%
  mutate(
    cluster = factor(as.character(cluster), levels = clusters_use),
    sample = as.character(sample)
  )

df_cell_long <- df_cell_wide %>%
  pivot_longer(
    cols = c(M1_score, M2_score),
    names_to = "signature",
    values_to = "score"
  ) %>%
  mutate(
    signature = recode(signature, M1_score = "M1", M2_score = "M2"),
    signature = factor(signature, levels = c("M1", "M2"))
  )

cell_summary <- df_cell_long %>%
  group_by(signature, cluster) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

cell_stats <- df_cell_long %>%
  group_by(signature) %>%
  summarise(
    n_cells_cluster0 = sum(cluster == clusters_use[[1]]),
    n_cells_cluster5 = sum(cluster == clusters_use[[2]]),
    mean_cluster0 = mean(score[cluster == clusters_use[[1]]], na.rm = TRUE),
    mean_cluster5 = mean(score[cluster == clusters_use[[2]]], na.rm = TRUE),
    median_cluster0 = median(score[cluster == clusters_use[[1]]], na.rm = TRUE),
    median_cluster5 = median(score[cluster == clusters_use[[2]]], na.rm = TRUE),
    test = "Wilcoxon rank-sum",
    p_value = safe_wilcox_p(score[cluster == clusters_use[[1]]], score[cluster == clusters_use[[2]]]),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    stars = p_stars(p_value)
  )

# ------------------------------------------------------------------------------
# SAMPLE-LEVEL DATA
# ------------------------------------------------------------------------------
df_sample <- df_cell_long %>%
  filter(!is.na(sample), sample != "") %>%
  group_by(sample, cluster, signature) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sample_short = recode(
      sample,
      "Hashtag-01-TotalSeqC" = "HTO1",
      "Hashtag-02-TotalSeqC" = "HTO2",
      "Hashtag-03-TotalSeqC" = "HTO3",
      "Hashtag-04-TotalSeqC" = "HTO4",
      "Hashtag-05-TotalSeqC" = "HTO5",
      "Hashtag-06-TotalSeqC" = "HTO6",
      "Hashtag-07-TotalSeqC" = "HTO7",
      "Hashtag-08-TotalSeqC" = "HTO8",
      .default = sample
    )
  )

sample_summary <- df_sample %>%
  group_by(signature, cluster) %>%
  summarise(
    n_samples = n(),
    mean_of_sample_means = mean(mean_score, na.rm = TRUE),
    median_of_sample_means = median(mean_score, na.rm = TRUE),
    sd_of_sample_means = sd(mean_score, na.rm = TRUE),
    mean_of_sample_medians = mean(median_score, na.rm = TRUE),
    median_of_sample_medians = median(median_score, na.rm = TRUE),
    .groups = "drop"
  )

paired_sample_table <- function(data, value_col) {
  data %>%
    select(sample, sample_short, cluster, signature, value = all_of(value_col)) %>%
    pivot_wider(names_from = cluster, values_from = value, names_prefix = "cluster_") %>%
    filter(!is.na(.data[[paste0("cluster_", clusters_use[[1]])]]), !is.na(.data[[paste0("cluster_", clusters_use[[2]])]])) %>%
    mutate(delta_0_minus_5 = .data[[paste0("cluster_", clusters_use[[1]])]] - .data[[paste0("cluster_", clusters_use[[2]])]])
}

sample_paired_mean <- paired_sample_table(df_sample, "mean_score")
sample_paired_median <- paired_sample_table(df_sample, "median_score")

sample_test <- function(data, signature_name, score_summary) {
  tmp <- data %>% filter(signature == signature_name)
  cluster0_col <- paste0("cluster_", clusters_use[[1]])
  cluster5_col <- paste0("cluster_", clusters_use[[2]])

  if (nrow(tmp) < 2) {
    return(data.frame(
      signature = signature_name,
      level = "Sample/HTO",
      score_summary = score_summary,
      n_paired_samples = nrow(tmp),
      cluster0 = NA_real_,
      cluster5 = NA_real_,
      mean_paired_delta_0_minus_5 = NA_real_,
      median_paired_delta_0_minus_5 = NA_real_,
      direction = NA_character_,
      test = "Paired Wilcoxon signed-rank",
      p_value = NA_real_
    ))
  }

  delta_mean <- mean(tmp$delta_0_minus_5, na.rm = TRUE)
  data.frame(
    signature = signature_name,
    level = "Sample/HTO",
    score_summary = score_summary,
    n_paired_samples = nrow(tmp),
    cluster0 = mean(tmp[[cluster0_col]], na.rm = TRUE),
    cluster5 = mean(tmp[[cluster5_col]], na.rm = TRUE),
    mean_paired_delta_0_minus_5 = delta_mean,
    median_paired_delta_0_minus_5 = median(tmp$delta_0_minus_5, na.rm = TRUE),
    direction = case_when(
      delta_mean > 0 ~ "Cluster 0 higher",
      delta_mean < 0 ~ "Cluster 5 higher",
      TRUE ~ "Equal"
    ),
    test = "Paired Wilcoxon signed-rank",
    p_value = safe_wilcox_p(tmp[[cluster0_col]], tmp[[cluster5_col]], paired = TRUE)
  )
}

sample_stats_mean <- bind_rows(
  sample_test(sample_paired_mean, "M1", "Mean score per cluster x HTO"),
  sample_test(sample_paired_mean, "M2", "Mean score per cluster x HTO")
) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

sample_stats_median <- bind_rows(
  sample_test(sample_paired_median, "M1", "Median score per cluster x HTO"),
  sample_test(sample_paired_median, "M2", "Median score per cluster x HTO")
) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

final_direction_summary <- bind_rows(
  cell_stats %>% transmute(signature, level = "Cell", direction = ifelse(mean_cluster0 > mean_cluster5, "Cluster 0 higher", "Cluster 5 higher"), p_value, p_adj_BH),
  sample_stats_mean %>% transmute(signature, level = "Sample mean", direction, p_value, p_adj_BH),
  sample_stats_median %>% transmute(signature, level = "Sample median", direction, p_value, p_adj_BH)
)

# ------------------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------------------
cols_fill <- c("M1_0" = "#ef4444", "M1_5" = "#f87171", "M2_0" = "#3b82f6", "M2_5" = "#93c5fd")

df_plot <- df_cell_long %>%
  mutate(sig_cluster = paste0(signature, "_", cluster))

med_df <- df_plot %>%
  group_by(signature, cluster) %>%
  summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    x = as.numeric(cluster),
    xstart = x - 0.12,
    xend = x + 0.12,
    sig_cluster = paste0(signature, "_", cluster)
  )

star_df <- df_plot %>%
  group_by(signature) %>%
  summarise(y = max(score, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  left_join(cell_stats %>% select(signature, stars), by = "signature") %>%
  mutate(x = 1.5)

p_violin <- ggplot(df_plot, aes(x = cluster, y = score, fill = sig_cluster, color = sig_cluster)) +
  geom_violin(width = 0.95, alpha = 0.65, linewidth = 0.6, trim = TRUE) +
  geom_segment(
    data = med_df,
    aes(x = xstart, xend = xend, y = median_score, yend = median_score, color = sig_cluster),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 1.1
  ) +
  geom_text(data = star_df, aes(x = x, y = y, label = stars), inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_wrap(~ signature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cols_fill) +
  scale_color_manual(values = cols_fill) +
  labs(title = "M1-like and M2-like signatures in clusters 0 and 5", x = "Cluster", y = "Module score") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"), panel.border = element_rect(fill = NA, linewidth = 0.6), plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(outdir, "M1_M2_cluster0_vs5_violin.pdf"), p_violin, width = 7.2, height = 4.5)
ggsave(file.path(outdir, "M1_M2_FINAL_violin_cluster0_RED_cluster5_BLUE.pdf"), p_violin, width = 7.2, height = 4.5)
ggsave(file.path(outdir, "M1_M2_FINAL_violin_cluster0_RED_cluster5_BLUE.png"), p_violin, width = 7.2, height = 4.5, dpi = 300)

p_ridge <- ggplot(df_plot, aes(x = score, y = cluster, fill = signature, color = signature)) +
  geom_density_ridges(alpha = 0.45, scale = 1.05, rel_min_height = 0.01, linewidth = 0.25, position = "identity") +
  facet_wrap(~ signature, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = c(M1 = "#ef4444", M2 = "#3b82f6")) +
  scale_color_manual(values = c(M1 = "#991b1b", M2 = "#1d4ed8")) +
  labs(title = "M1-like and M2-like signatures: clusters 0 and 5", x = "Module score", y = "Cluster") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "M1_M2_clusters_0_5_overlay_ridgeplot.pdf"), p_ridge, width = 8, height = 5)
ggsave(file.path(outdir, "M1_M2_clusters_0_5_overlay_ridgeplot.png"), p_ridge, width = 8, height = 5, dpi = 300)

p_sample <- ggplot(df_sample, aes(x = cluster, y = mean_score, group = sample_short, color = sample_short)) +
  geom_line(alpha = 0.55, linewidth = 0.4) +
  geom_point(size = 2) +
  facet_wrap(~ signature, nrow = 1, scales = "free_y") +
  labs(title = "M1-like and M2-like signatures: paired sample-level comparison", x = "Cluster", y = "Mean module score per cluster x sample") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "M1_M2_clusters_0_5_sample_level_paired.pdf"), p_sample, width = 8, height = 4.8)
ggsave(file.path(outdir, "M1_M2_clusters_0_5_sample_level_paired.png"), p_sample, width = 8, height = 4.8, dpi = 300)

p_box <- ggplot(df_sample, aes(x = cluster, y = mean_score, fill = cluster, color = cluster)) +
  geom_boxplot(alpha = 0.45, outlier.shape = NA, width = 0.55) +
  geom_point(aes(group = sample_short), position = position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.8) +
  facet_wrap(~ signature, nrow = 1, scales = "free_y") +
  labs(title = "M1-like and M2-like signatures at sample level", x = "Cluster", y = "Mean module score per cluster x sample") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "M1_M2_clusters_0_5_sample_level_boxplot.pdf"), p_box, width = 7.2, height = 4.8)
ggsave(file.path(outdir, "M1_M2_clusters_0_5_sample_level_boxplot.png"), p_box, width = 7.2, height = 4.8, dpi = 300)

if ("umap" %in% Reductions(obj)) {
  df_umap <- FetchData(obj, vars = c(m1_score$column, m2_score$column, cluster_col)) %>%
    mutate(
      M1_M2_delta = .data[[m1_score$column]] - .data[[m2_score$column]],
      cluster = as.character(.data[[cluster_col]])
    )
  emb <- Embeddings(obj, reduction = "umap") %>% as.data.frame()
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  df_umap <- bind_cols(emb, df_umap)

  p_umap <- ggplot(df_umap, aes(UMAP_1, UMAP_2)) +
    geom_point(data = df_umap %>% filter(!cluster %in% selected_umap_clusters), color = "grey82", size = 0.12) +
    geom_point(data = df_umap %>% filter(cluster %in% selected_umap_clusters), aes(color = M1_M2_delta), size = 0.18) +
    scale_color_gradient2(low = "#2563eb", mid = "grey90", high = "#dc2626", midpoint = 0, name = "M1 - M2") +
    labs(title = "M1 / M2 signature", subtitle = paste("Colored clusters:", paste(selected_umap_clusters, collapse = ", "))) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path(outdir, "M1_M2_selected_clusters.pdf"), p_umap, width = 6.5, height = 5.5)
  ggsave(file.path(outdir, "M1_M2_selected_clusters.svg"), p_umap, width = 6.5, height = 5.5)
  ggsave(file.path(outdir, "M1_M2_selected_clusters.jpg"), p_umap, width = 6.5, height = 5.5, dpi = 300)
}

# ------------------------------------------------------------------------------
# EXPORTS
# ------------------------------------------------------------------------------
write.csv(final_direction_summary, file.path(outdir, "M1_M2_FINAL_direction_summary.csv"), row.names = FALSE)
write.csv(cell_stats, file.path(outdir, "M1_M2_cell_level_statistics_cluster0_vs_5.csv"), row.names = FALSE)
write.csv(cell_summary, file.path(outdir, "M1_M2_cell_level_summary.csv"), row.names = FALSE)
write.csv(df_sample, file.path(outdir, "M1_M2_cluster_by_HTO_scores.csv"), row.names = FALSE)
write.csv(sample_stats_mean, file.path(outdir, "M1_M2_sample_level_primary_statistics.csv"), row.names = FALSE)
write.csv(sample_stats_median, file.path(outdir, "M1_M2_sample_level_median_robustness_statistics.csv"), row.names = FALSE)

workbook <- list(
  FINAL_direction = final_direction_summary,
  Gene_presence = gene_check,
  Cell_stats = cell_stats,
  Cell_summary = cell_summary,
  Sample_primary_stats = sample_stats_mean,
  Sample_median_check = sample_stats_median,
  Sample_summary = sample_summary,
  Paired_sample_means = sample_paired_mean,
  Paired_sample_medians = sample_paired_median,
  Cluster_HTO_scores = df_sample,
  Cell_level_scores = df_cell_long
)

openxlsx::write.xlsx(workbook, file.path(outdir, "M1_M2_FINAL_all_statistics.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(list(cell_level_stats = cell_stats, hashID_level_stats = sample_stats_mean), file.path(outdir, "Stats_cluster0_vs5.xlsx"), overwrite = TRUE)

capture.output(sessionInfo(), file = file.path(outdir, "M1_M2_sessionInfo.txt"))

if (save_scored_rds) {
  saveRDS(obj, file.path(outdir, "seurat_final_clean4_M1_M2_scored.rds"))
}

message("Done. Outputs written to: ", outdir)
