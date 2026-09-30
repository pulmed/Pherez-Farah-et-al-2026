#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/qc/biotin_delta_threshold_sensitivity.R
# Original files: 20260222 NOISE FLOOR DELTA.R; 20260222 BOOTSTRAP DELTA.R; 20260222 COHEN DELTA.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Evaluate the sample-level Biotin delta cutoff with noise-floor, bootstrap, and effect-size summaries.
# Inputs: Canonical analysis-ready Seurat object with ADT data plus final_clusters, condition, and hash.ID metadata.
# Outputs: Sample-level Biotin delta sensitivity tables, workbooks, and plots.
# Assay/layer input: ADT data.
# Dependencies: Seurat, dplyr, tidyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - The bootstrap/noise-floor unit is the hashtag/sample, not the individual cell.
# - The 0.1 delta cutoff is treated as a practical effect-size guide.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

set.seed(as.integer(Sys.getenv("SEED", unset = "1")))

rds_file <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_dir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "05_biotin_analysis", "delta_threshold_sensitivity")

cluster_col <- "final_clusters"
condition_col <- "condition"
sample_col <- "hash.ID"
adt_assay <- "ADT"
biotin_feature <- "Biotin-TotalSeqC"

treated_label <- "treated"
untreated_label <- "untreated"
clusters_of_interest <- trimws(strsplit(
  Sys.getenv("CLUSTERS_OF_INTEREST", unset = "10,0,1,2,3,4,5,7,8,9,13,14"),
  ",",
  fixed = TRUE
)[[1]])

min_cells_per_sample_cluster <- as.integer(Sys.getenv("MIN_CELLS_PER_SAMPLE_CLUSTER", unset = "20"))
n_iter <- as.integer(Sys.getenv("N_ITER", unset = "10000"))
n_boot <- as.integer(Sys.getenv("N_BOOT", unset = "20000"))
delta_cut <- as.numeric(Sys.getenv("DELTA_CUT", unset = "0.1"))
conf_level <- as.numeric(Sys.getenv("CONF_LEVEL", unset = "0.95"))

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(rds_file)) {
  message("SMOKE_TEST=1: Seurat object not found; skipping Biotin delta sensitivity.")
  quit(save = "no", status = 0)
}

if (!file.exists(rds_file)) stop("Cannot find Seurat object: ", rds_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(rds_file)
missing_meta <- setdiff(c(cluster_col, condition_col, sample_col), colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
if (!adt_assay %in% Assays(obj)) stop("Missing assay: ", adt_assay)
if (!biotin_feature %in% rownames(obj[[adt_assay]])) stop("Missing ADT feature: ", biotin_feature)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: Biotin feature present: ", biotin_feature)
  message("SMOKE_TEST=1: required metadata present: ", paste(c(cluster_col, condition_col, sample_col), collapse = ", "))
  message("SMOKE_TEST=1: skipping full Biotin delta sensitivity analysis.")
  quit(save = "no", status = 0)
}

DefaultAssay(obj) <- adt_assay
obj <- NormalizeData(obj, assay = adt_assay, normalization.method = "CLR", margin = 2, verbose = FALSE)

condition <- tolower(trimws(as.character(obj@meta.data[[condition_col]])))
condition <- dplyr::recode(
  condition,
  "control" = untreated_label,
  "ctrl" = untreated_label,
  "vehicle" = untreated_label,
  "mock" = untreated_label,
  "untreated" = untreated_label,
  "treated" = treated_label,
  .default = condition
)
obj@meta.data[[condition_col]] <- factor(condition, levels = c(untreated_label, treated_label))

metadata <- obj@meta.data %>%
  mutate(
    cell = rownames(obj@meta.data),
    cluster = as.character(.data[[cluster_col]]),
    condition = as.character(.data[[condition_col]]),
    sample_id = as.character(.data[[sample_col]])
  )

biotin_values <- as.numeric(GetAssayData(obj, assay = adt_assay, slot = "data")[biotin_feature, metadata$cell])

df_cell <- metadata %>%
  mutate(biotin = biotin_values) %>%
  filter(condition %in% c(untreated_label, treated_label), !is.na(cluster), !is.na(sample_id))

df_sample <- df_cell %>%
  group_by(cluster, condition, sample_id) %>%
  summarise(
    n_cells = n(),
    sample_median = median(biotin),
    sample_mean = mean(biotin),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample_cluster)

if (length(clusters_of_interest) > 0 && !all(is.na(clusters_of_interest)) && any(nzchar(clusters_of_interest))) {
  df_sample <- df_sample %>% filter(cluster %in% clusters_of_interest)
}

observed_delta <- df_sample %>%
  group_by(cluster) %>%
  summarise(
    n_untreated = sum(condition == untreated_label),
    n_treated = sum(condition == treated_label),
    median_untreated = ifelse(n_untreated > 0, median(sample_median[condition == untreated_label]), NA_real_),
    median_treated = ifelse(n_treated > 0, median(sample_median[condition == treated_label]), NA_real_),
    delta_median = median_treated - median_untreated,
    .groups = "drop"
  )

split_null_delta <- function(x, n_iter) {
  n <- length(x)
  if (n < 4) return(rep(NA_real_, n_iter))

  half <- floor(n / 2)
  replicate(n_iter, {
    idx <- sample.int(n, size = half, replace = FALSE)
    median(x[idx]) - median(x[-idx])
  })
}

null_groups <- df_sample %>%
  group_by(cluster, condition) %>%
  summarise(values = list(sample_median), n_samples = n(), .groups = "drop")

null_groups$null_delta <- lapply(null_groups$values, split_null_delta, n_iter = n_iter)

null_long <- null_groups %>%
  select(cluster, condition, n_samples, null_delta) %>%
  tidyr::unnest_longer(null_delta, values_to = "delta_null")

noise_floor <- null_long %>%
  group_by(cluster, condition) %>%
  summarise(
    n_samples = first(n_samples),
    q95_abs = quantile(abs(delta_null), 0.95, na.rm = TRUE),
    q99_abs = quantile(abs(delta_null), 0.99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cluster) %>%
  summarise(
    q95_abs_max = max(q95_abs, na.rm = TRUE),
    q99_abs_max = max(q99_abs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(observed_delta, by = "cluster") %>%
  mutate(
    passes_delta_cut = !is.na(delta_median) & delta_median > delta_cut,
    exceeds_noise_q95 = !is.na(delta_median) & abs(delta_median) > q95_abs_max,
    exceeds_noise_q99 = !is.na(delta_median) & abs(delta_median) > q99_abs_max
  )

boot_delta_median <- function(x_untreated, x_treated, n_boot) {
  n_untreated <- length(x_untreated)
  n_treated <- length(x_treated)
  if (n_untreated < 2 || n_treated < 2) return(rep(NA_real_, n_boot))

  replicate(n_boot, {
    median(sample(x_treated, size = n_treated, replace = TRUE)) -
      median(sample(x_untreated, size = n_untreated, replace = TRUE))
  })
}

alpha <- 1 - conf_level
bootstrap_delta <- lapply(sort(unique(df_sample$cluster)), function(cluster_id) {
  dat <- df_sample %>% filter(cluster == cluster_id)
  x_untreated <- dat$sample_median[dat$condition == untreated_label]
  x_treated <- dat$sample_median[dat$condition == treated_label]
  boot <- boot_delta_median(x_untreated, x_treated, n_boot = n_boot)

  data.frame(
    cluster = cluster_id,
    n_untreated = length(x_untreated),
    n_treated = length(x_treated),
    delta_median = ifelse(
      length(x_untreated) > 0 && length(x_treated) > 0,
      median(x_treated) - median(x_untreated),
      NA_real_
    ),
    ci_low = as.numeric(quantile(boot, probs = alpha / 2, na.rm = TRUE)),
    ci_high = as.numeric(quantile(boot, probs = 1 - alpha / 2, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
})
bootstrap_delta <- bind_rows(bootstrap_delta) %>%
  mutate(
    passes_delta_cut = !is.na(delta_median) & delta_median > delta_cut,
    ci_excludes_zero = !is.na(ci_low) & !is.na(ci_high) & (ci_low > 0 | ci_high < 0),
    ci_above_delta_cut = !is.na(ci_low) & ci_low > delta_cut
  )

effect_size <- function(dat) {
  untreated <- dat$sample_median[dat$condition == untreated_label]
  treated <- dat$sample_median[dat$condition == treated_label]
  n_untreated <- length(untreated)
  n_treated <- length(treated)

  pooled_sd <- sqrt(((n_untreated - 1) * stats::var(untreated) + (n_treated - 1) * stats::var(treated)) /
    (n_untreated + n_treated - 2))
  cohen_d <- ifelse(is.finite(pooled_sd) && pooled_sd > 0, (mean(treated) - mean(untreated)) / pooled_sd, NA_real_)
  correction <- 1 - (3 / (4 * (n_untreated + n_treated) - 9))

  data.frame(
    n_untreated = n_untreated,
    n_treated = n_treated,
    median_untreated = median(untreated),
    median_treated = median(treated),
    delta_median = median(treated) - median(untreated),
    cohen_d = cohen_d,
    hedges_g = cohen_d * correction,
    p_wilcox = tryCatch(stats::wilcox.test(sample_median ~ condition, data = dat, exact = FALSE)$p.value, error = function(e) NA_real_)
  )
}

effect_stats <- df_sample %>%
  group_by(cluster) %>%
  filter(
    sum(condition == untreated_label) >= 2,
    sum(condition == treated_label) >= 2
  ) %>%
  group_modify(~ effect_size(.x)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH"))

write.csv(df_sample, file.path(output_dir, "biotin_sample_medians.csv"), row.names = FALSE)
write.csv(observed_delta, file.path(output_dir, "observed_delta_by_cluster.csv"), row.names = FALSE)
write.csv(noise_floor, file.path(output_dir, "noise_floor_vs_observed_delta.csv"), row.names = FALSE)
write.csv(bootstrap_delta, file.path(output_dir, "bootstrap_delta_ci_by_cluster.csv"), row.names = FALSE)
write.csv(effect_stats, file.path(output_dir, "effect_size_by_cluster.csv"), row.names = FALSE)

write.xlsx(
  list(
    sample_medians = df_sample,
    observed_delta = observed_delta,
    noise_floor = noise_floor,
    bootstrap_delta = bootstrap_delta,
    effect_size = effect_stats
  ),
  file.path(output_dir, "biotin_delta_threshold_sensitivity.xlsx"),
  overwrite = TRUE
)

plot_delta <- bootstrap_delta %>%
  arrange(delta_median) %>%
  mutate(cluster = factor(cluster, levels = cluster))

p_ci <- ggplot(plot_delta, aes(x = delta_median, y = cluster)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = delta_cut, linetype = "dotted", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("Delta median, treated - untreated; dotted line = ", delta_cut),
    y = "Cluster"
  )

ggsave(file.path(output_dir, "plots", "bootstrap_delta_ci_by_cluster.pdf"), p_ci, width = 8, height = 6)

p_noise <- noise_floor %>%
  arrange(delta_median) %>%
  mutate(cluster = factor(cluster, levels = cluster)) %>%
  ggplot(aes(x = abs(delta_median), y = cluster)) +
  geom_vline(xintercept = delta_cut, linetype = "dotted", linewidth = 0.45) +
  geom_point(size = 2) +
  geom_point(aes(x = q95_abs_max), shape = 1, size = 2) +
  theme_classic(base_size = 12) +
  labs(
    x = "Absolute observed delta (filled) and 95% null split delta (open)",
    y = "Cluster"
  )

ggsave(file.path(output_dir, "plots", "noise_floor_vs_observed_delta.pdf"), p_noise, width = 8, height = 6)

writeLines(
  c(
    "Biotin delta threshold sensitivity settings",
    paste0("Seurat object: ", rds_file),
    paste0("ADT feature: ", biotin_feature),
    paste0("Minimum cells per sample-cluster: ", min_cells_per_sample_cluster),
    paste0("Noise-floor iterations: ", n_iter),
    paste0("Bootstrap iterations: ", n_boot),
    paste0("Delta cutoff: ", delta_cut),
    "Unit of inference: sample/hashtag-level Biotin medians."
  ),
  file.path(output_dir, "settings_used.txt")
)

message("Done. Biotin delta sensitivity outputs written to: ", output_dir)
