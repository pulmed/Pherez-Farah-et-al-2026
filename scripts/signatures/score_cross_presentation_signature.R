#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/score_cross_presentation_signature.R
# Original file: 20260823 CROOS PRES FINAL.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Score and summarize the final cross-presentation signature across mouse clusters.
# Inputs: Canonical analysis-ready Seurat object with RNA data, cluster, condition, and sample metadata.
# Outputs: Gene checks, ridgeplots, condition statistics, and sample-level summaries.
# Assay/layer input: RNA data.
# Dependencies: Seurat, SeuratObject, dplyr, tidyr, ggplot2, ggridges, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Uses the 26-gene 2026-08-23 signature with Wdfy4, Sec22b, and Lnpep additions.
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
outdir <- make_output_dir(file.path(output_root, "17_cross_presentation_signature"))

assay_use <- Sys.getenv("ASSAY", unset = "RNA")
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "final_clusters")
condition_col <- Sys.getenv("CONDITION_COL", unset = "condition")
sample_col_candidates <- parse_csv_env("SAMPLE_COL_CANDIDATES", c("HTO_classification", "hash.ID"))
clusters_use <- parse_csv_env("CLUSTERS_USE", c("0", "2", "3", "4", "5", "14", "13", "1", "7", "8", "9"))
condition_levels <- parse_csv_env("CONDITION_LEVELS", c("untreated", "treated"))
condition_labels <- parse_csv_env("CONDITION_LABELS", c("Untreated", "Treated"))
score_seed <- as.integer(Sys.getenv("MODULE_SCORE_SEED", unset = "1234"))
save_scored_rds <- identical(tolower(Sys.getenv("SAVE_SCORED_RDS", unset = "false")), "true")

cross_presentation_genes <- c(
  "H2-K1", "H2-D1", "B2m",
  "Tap1", "Tap2", "Tapbp",
  "Calr", "Canx", "Pdia3", "Erap1",
  "Psmb8", "Psmb9", "Psmb10", "Psme1", "Psme2",
  "Sec61a1", "Sec61b", "Sec61g",
  "Cybb", "Rab27a", "Rac2", "Vamp8", "Stx4",
  "Wdfy4", "Sec22b", "Lnpep"
)

# ------------------------------------------------------------------------------
# LOAD AND SCORE
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_assays(obj, assay_use)
require_metadata(obj, c(cluster_col, condition_col))
sample_col <- first_existing_metadata(obj, sample_col_candidates)

DefaultAssay(obj) <- assay_use
obj <- coerce_metadata_character(obj, c(cluster_col, condition_col, sample_col))
obj <- ensure_rna_data_layer(obj, assay = assay_use)

features <- rownames(obj[[assay_use]])
present_genes <- present_signature_genes(cross_presentation_genes, features, "CrossPresentation")
gene_check <- gene_presence_table(cross_presentation_genes, features, "CrossPresentation")

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: cross-presentation genes present: ", length(present_genes), "/", length(cross_presentation_genes))
  message("SMOKE_TEST=1: sample column selected: ", sample_col)
  quit(save = "no", status = 0)
}

write.csv(
  gene_check,
  file.path(outdir, "CrossPresentation_signature_gene_presence.csv"),
  row.names = FALSE
)

score <- add_module_score(obj, present_genes, "CrossPresentationSig_v2", assay = assay_use, seed = score_seed)
obj <- score$object
score_col <- score$column

# ------------------------------------------------------------------------------
# CELL-LEVEL CLUSTER SUMMARY
# ------------------------------------------------------------------------------
plot_df <- FetchData(obj, vars = c(score_col, cluster_col, condition_col, sample_col)) %>%
  rename(
    score = all_of(score_col),
    cluster = all_of(cluster_col),
    condition = all_of(condition_col),
    sample = all_of(sample_col)
  ) %>%
  mutate(
    cluster = as.character(cluster),
    condition = as.character(condition),
    sample = as.character(sample)
  ) %>%
  filter(cluster %in% clusters_use)

cluster_order <- plot_df %>%
  group_by(cluster) %>%
  summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_score)) %>%
  pull(cluster)

plot_df <- plot_df %>%
  mutate(cluster = factor(cluster, levels = rev(cluster_order)))

cluster_summary <- plot_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_score))

write.csv(
  cluster_summary,
  file.path(outdir, "CrossPresentation_FINAL_cluster_summary.csv"),
  row.names = FALSE
)

p_cluster <- ggplot(plot_df, aes(x = score, y = cluster)) +
  geom_density_ridges(fill = "#4db6ac", color = "#00695c", alpha = 0.55, scale = 1.1, rel_min_height = 0.01, linewidth = 0.25) +
  labs(
    title = "Cross-presentation signature across clusters",
    subtitle = paste("Clusters ordered high to low median:", paste(cluster_order, collapse = ", ")),
    x = "Module score",
    y = "Cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "CrossPresentation_FINAL_ridgeplot_high_to_low.png"), p_cluster, width = 10, height = 7, dpi = 300)
ggsave(file.path(outdir, "CrossPresentation_FINAL_ridgeplot_high_to_low.pdf"), p_cluster, width = 10, height = 7)

# ------------------------------------------------------------------------------
# CONDITION SUMMARY
# ------------------------------------------------------------------------------
condition_df <- plot_df %>%
  filter(condition %in% condition_levels) %>%
  mutate(condition = factor(condition, levels = condition_levels, labels = condition_labels))

condition_stats <- condition_df %>%
  group_by(cluster) %>%
  summarise(
    n_untreated = sum(condition == condition_labels[[1]]),
    n_treated = sum(condition == condition_labels[[2]]),
    mean_untreated = mean(score[condition == condition_labels[[1]]], na.rm = TRUE),
    mean_treated = mean(score[condition == condition_labels[[2]]], na.rm = TRUE),
    median_untreated = median(score[condition == condition_labels[[1]]], na.rm = TRUE),
    median_treated = median(score[condition == condition_labels[[2]]], na.rm = TRUE),
    delta_mean_treated_minus_untreated = mean_treated - mean_untreated,
    delta_median_treated_minus_untreated = median_treated - median_untreated,
    wilcox_p = safe_wilcox_p(score[condition == condition_labels[[1]]], score[condition == condition_labels[[2]]]),
    .groups = "drop"
  ) %>%
  mutate(wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")) %>%
  arrange(match(as.character(cluster), rev(levels(condition_df$cluster))))

openxlsx::write.xlsx(
  list(
    cluster_condition_stats = condition_stats,
    cell_level_scores = condition_df,
    signature_gene_check = gene_check
  ),
  file.path(outdir, "CrossPresentation_FINAL_condition_stats_by_cluster.xlsx"),
  overwrite = TRUE
)

p_condition <- ggplot(condition_df, aes(x = score, y = cluster, fill = condition, color = condition)) +
  geom_density_ridges(alpha = 0.45, scale = 1.15, rel_min_height = 0.01, linewidth = 0.25, position = "identity") +
  scale_fill_manual(values = setNames(c("#A6CEE3", "#F4A6A6"), condition_labels)) +
  scale_color_manual(values = setNames(c("#6BAED6", "#E57373"), condition_labels)) +
  labs(title = "Cross-presentation signature across clusters by condition", x = "Module score", y = "Cluster", fill = "Condition", color = "Condition") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(outdir, "CrossPresentation_FINAL_condition_ridgeplot.png"), p_condition, width = 10, height = 7, dpi = 300)
ggsave(file.path(outdir, "CrossPresentation_FINAL_condition_ridgeplot.pdf"), p_condition, width = 10, height = 7)

# ------------------------------------------------------------------------------
# SAMPLE-LEVEL SUMMARY
# ------------------------------------------------------------------------------
df_rep <- plot_df %>%
  filter(!is.na(sample), sample != "") %>%
  group_by(cluster, sample) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  df_rep,
  file.path(outdir, "CrossPresentation_FINAL_cluster_by_hashtag_replicates.csv"),
  row.names = FALSE
)

replicate_cluster_summary <- df_rep %>%
  group_by(cluster) %>%
  summarise(
    n_samples = n(),
    mean_of_sample_means = mean(mean_score, na.rm = TRUE),
    median_of_sample_means = median(mean_score, na.rm = TRUE),
    sd_of_sample_means = sd(mean_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_of_sample_means))

write.csv(
  replicate_cluster_summary,
  file.path(outdir, "CrossPresentation_FINAL_sample_level_cluster_summary.csv"),
  row.names = FALSE
)

p_rep <- ggplot(df_rep, aes(x = cluster, y = mean_score, fill = cluster)) +
  geom_boxplot(alpha = 0.45, outlier.shape = NA, width = 0.55) +
  geom_point(aes(color = sample), position = position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.85) +
  labs(title = "Cross-presentation signature across clusters", subtitle = "Each point is one hashtag/sample replicate within a cluster", x = "Cluster", y = "Mean module score per cluster x sample") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "CrossPresentation_FINAL_cluster_hashtag_replicates_boxplot.png"), p_rep, width = 10, height = 6, dpi = 300)
ggsave(file.path(outdir, "CrossPresentation_FINAL_cluster_hashtag_replicates_boxplot.pdf"), p_rep, width = 10, height = 6)

if (save_scored_rds) {
  saveRDS(obj, file.path(outdir, "seurat_final_clean4_CrossPresentation_FINAL_scored.rds"))
}

message("Done. Outputs written to: ", outdir)
