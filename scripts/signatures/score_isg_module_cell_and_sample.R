#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/score_isg_module_cell_and_sample.R
# Original file: 20260218 ISG CELL AND SAMPLE.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Score ISG modules at cell and sample level.
# Inputs: Canonical analysis-ready Seurat object with RNA data and sample metadata.
# Outputs: ISG cell-level and sample-level workbooks and plots.
# Assay/layer input: RNA data.
# Dependencies: Seurat, dplyr, ggplot2, openxlsx, pheatmap.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Writes outputs under output/12_20260218_ISG_MODULE by default.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(stringr)
  library(rlang)
})

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "12_20260218_ISG_MODULE")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cluster_col <- "final_clusters"
hash_col <- "hash.ID"
exclude_clusters <- c("17")

# Panel cluster orders
clusters_A <- c("1","7","8","9","13","12","10")   # 7 clusters
clusters_B <- c("0","2","3","4","5","14")         # 6 clusters (12/10 removed; 17 ignored)

# colors: control first (grey), treated second (blue)
cond_levels <- c("untreated","treated")
cond_colors <- c(untreated = "grey80", treated = "#7C9FE6")

# signature genes
sig_genes <- c(
  "Irgm2","Fcgr3","Il10ra","Pik3ap1","Oas3","Ifi47","Entpd1","Tnfaip2","Rtp4","Gna13",
  "Irgm1","Fcgr4","Stat2","Sdcbp","Gm4951","Ifi213","Samd9l","Igtp","Ifi207","Xaf1",
  "Socs3","Mndal","Pde7b","Irf7","Gm12185","Slamf8","Atp11b","Samhd1","Fbxl5","Tgtp2",
  "Slfn1","Cd86","Ifi27l2a","Gm20663","Zbp1","Slfn8","Creb5","Irak2","Ifi206","Pnp",
  "Iigp1","Ly6i","Batf2","Arid5b","Bach1","Ifi204","Oasl2","Ifi203","Parp14","Gbp3",
  "Slfn2","Ptprc","Cyrib","Tgfbi","Gpr141","Ifit2","Cd300lf","Gbp7","Rnf213","Fcgr1",
  "Fgl2","Atp6v0c","Ifi205","Gbp4","Cebpb","Lilrb4b","Tgm2","Pkm","Gbp2","Ifi211",
  "Lilrb4a","Cxcl10","Ctss","Ly6a"
)

# AddModuleScore naming
score_base <- "ISGsig"
score_col <- paste0(score_base, "1")

# ------------------------------------------------------------------------------
# YOUR BASE STYLE (A) - THINNER
# ------------------------------------------------------------------------------
# (These were the only numbers changed to make BOTH panels thin like your preferred A)
nudge_A <- 0.04
alpha_vals <- c(untreated = 0.50, treated = 0.95)
violin_width_A <- 0.70
adjust_k <- 1.2

median_half_A <- 0.07
median_lwd_u <- 0.9
median_lwd_t <- 1.1
median_col_u <- "grey30"
median_col_t <- "#1E4F9E"
median_lty <- "11"

# ------------------------------------------------------------------------------
# MAKE B LOOK IDENTICAL TO A
# ------------------------------------------------------------------------------
# B has fewer clusters -> if saved at same canvas width, violins look fatter.
# Scale B geometry by nB/nA to match A's perceived thickness.
scale_B <- length(clusters_B) / length(clusters_A)   # 6/7

nudge_B <- nudge_A        * scale_B
violin_width_B <- violin_width_A * scale_B
median_half_B <- median_half_A  * scale_B

# ------------------------------------------------------------------------------
# PRECHECKS
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
stopifnot(exists("seurat_final_clean4"))
stopifnot(all(c(cluster_col, hash_col) %in% colnames(seurat_final_clean4[[]])))

DefaultAssay(seurat_final_clean4) <- "RNA"
if (ncol(GetAssayData(seurat_final_clean4, slot = "data")) == 0) {
  seurat_final_clean4 <- NormalizeData(seurat_final_clean4, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# CONDITION FROM HASH.ID
# ------------------------------------------------------------------------------
hash_raw <- as.character(seurat_final_clean4[[hash_col]][,1])
hash_num <- stringr::str_match(hash_raw, "Hashtag-(\\d+)-")[,2]
hash_int <- suppressWarnings(as.integer(hash_num))

seurat_final_clean4$condition <- ifelse(hash_int %in% 1:4, "treated",
                                        ifelse(hash_int %in% 5:8, "untreated", NA))
seurat_final_clean4$condition <- factor(seurat_final_clean4$condition, levels = cond_levels)

message("Condition counts (after mapping):")
print(table(seurat_final_clean4$condition, useNA="ifany"))

# ------------------------------------------------------------------------------
# MODULE SCORE (compute or reuse)
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(seurat_final_clean4))
missing <- setdiff(sig_genes, present)
message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing genes ignored: ", paste(missing, collapse = ", "))

if (!score_col %in% colnames(seurat_final_clean4[[]])) {
  seurat_final_clean4 <- AddModuleScore(
    seurat_final_clean4,
    features = list(present),
    name = score_base,
    verbose = FALSE
  )
}
message("Using score column: ", score_col)

# ------------------------------------------------------------------------------
# CELL-LEVEL DF
# ------------------------------------------------------------------------------
df_cell <- FetchData(seurat_final_clean4, vars = c(score_col, cluster_col, "condition", hash_col)) %>%
  dplyr::rename(
    score = !!sym(score_col),
    cluster = !!sym(cluster_col),
    sample = !!sym(hash_col)
  ) %>%
  mutate(
    cluster = as.character(cluster),
    sample = as.character(sample),
    condition = as.character(condition)
  ) %>%
  filter(condition %in% cond_levels) %>%
  filter(!cluster %in% exclude_clusters)

df_cell$condition <- factor(df_cell$condition, levels = cond_levels)

# ------------------------------------------------------------------------------
# SAMPLE-LEVEL DF (median per sample x cluster)
# ------------------------------------------------------------------------------
df_sample <- df_cell %>%
  group_by(sample, condition, cluster) %>%
  summarise(
    n_cells = n(),
    score = median(score, na.rm = TRUE),
    .groups = "drop"
  )
df_sample$condition <- factor(df_sample$condition, levels = cond_levels)

# ------------------------------------------------------------------------------
# OVERLAY VIOLIN FUNCTION (parametrized)
# ------------------------------------------------------------------------------
plot_overlay <- function(d_in, clusters_vec, title_text, subtitle_text, ylab_text,
                         nudge_val, violin_width, median_half_width) {

  d <- d_in %>% filter(cluster %in% clusters_vec)
  stopifnot(nrow(d) > 0)

  d$cluster <- factor(d$cluster, levels = clusters_vec)

  # robust clipping (keeps y scale comparable + removes extreme tails)
  lo <- as.numeric(quantile(d$score, 0.05, na.rm = TRUE))
  hi <- as.numeric(quantile(d$score, 0.95, na.rm = TRUE))
  d$score_clip <- pmin(pmax(d$score, lo), hi)

  # numeric x for manual dodge
  d$x_base <- as.numeric(d$cluster)
  d$x_ut <- d$x_base - nudge_val
  d$x_tr <- d$x_base + nudge_val

  # median "stubs"
  med_tab <- aggregate(score_clip ~ cluster + condition, d, median, na.rm = TRUE)
  med_tab$x <- ifelse(med_tab$condition == "untreated",
                      as.numeric(med_tab$cluster) - nudge_val,
                      as.numeric(med_tab$cluster) + nudge_val)
  med_tab$x_min <- med_tab$x - median_half_width
  med_tab$x_max <- med_tab$x + median_half_width
  med_tab$y <- med_tab$score_clip

  ggplot() +
    # ------------------------------------------------------------------------------
    # untreated violin
    # ------------------------------------------------------------------------------
  geom_violin(
    data = d[d$condition == "untreated", ],
    aes(x = x_ut, y = score_clip, group = cluster),
    fill = cond_colors["untreated"], color = "grey50",
    width = violin_width, adjust = adjust_k, alpha = alpha_vals["untreated"],
    linewidth = 0.5, trim = FALSE,
    scale = "width"   # FIX 1
  ) +
    # ------------------------------------------------------------------------------
    # treated violin
    # ------------------------------------------------------------------------------
  geom_violin(
    data = d[d$condition == "treated", ],
    aes(x = x_tr, y = score_clip, group = cluster),
    fill = cond_colors["treated"], color = "#335C99",
    width = violin_width, adjust = adjust_k, alpha = alpha_vals["treated"],
    linewidth = 0.7, trim = FALSE,
    scale = "width"   # FIX 1
  ) +
    # ------------------------------------------------------------------------------
    # medians (untreated)
    # ------------------------------------------------------------------------------
  geom_segment(
    data = subset(med_tab, condition == "untreated"),
    aes(x = x_min, xend = x_max, y = y, yend = y),
    linewidth = median_lwd_u, linetype = median_lty,
    color = median_col_u, lineend = "round"
  ) +
    # ------------------------------------------------------------------------------
    # medians (treated)
    # ------------------------------------------------------------------------------
  geom_segment(
    data = subset(med_tab, condition == "treated"),
    aes(x = x_min, xend = x_max, y = y, yend = y),
    linewidth = median_lwd_t, linetype = median_lty,
    color = median_col_t, lineend = "round"
  ) +
    scale_x_continuous(
      breaks = seq_along(clusters_vec),
      labels = clusters_vec,
      expand = c(0, 0)
    ) +
    labs(x = "Cluster", y = ylab_text, title = title_text, subtitle = subtitle_text) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.ticks.length.x = grid::unit(2, "pt"),
      plot.title = element_text(face = "bold")
    )
}

# ------------------------------------------------------------------------------
# STATS TABLES
# ------------------------------------------------------------------------------
safe_wilcox <- function(d) {
  if (!all(cond_levels %in% d$condition)) return(NA_real_)
  out <- try(wilcox.test(score ~ condition, data = d, exact = FALSE)$p.value, silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

cell_stats <- df_cell %>%
  group_by(cluster) %>%
  summarise(
    n_untreated = sum(condition=="untreated"),
    n_treated = sum(condition=="treated"),
    median_untreated = median(score[condition=="untreated"], na.rm=TRUE),
    median_treated = median(score[condition=="treated"],   na.rm=TRUE),
    mean_untreated = mean(score[condition=="untreated"], na.rm=TRUE),
    mean_treated = mean(score[condition=="treated"],   na.rm=TRUE),
    delta_median = median_treated - median_untreated,
    delta_mean = mean_treated - mean_untreated,
    p_value = safe_wilcox(cur_data()),
    .groups="drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method="BH")) %>%
  arrange(suppressWarnings(as.numeric(cluster)), cluster)

sample_stats <- df_sample %>%
  group_by(cluster) %>%
  summarise(
    n_samples_untreated = n_distinct(sample[condition=="untreated"]),
    n_samples_treated = n_distinct(sample[condition=="treated"]),
    median_untreated = median(score[condition=="untreated"], na.rm=TRUE),
    median_treated = median(score[condition=="treated"],   na.rm=TRUE),
    mean_untreated = mean(score[condition=="untreated"], na.rm=TRUE),
    mean_treated = mean(score[condition=="treated"],   na.rm=TRUE),
    delta_median = median_treated - median_untreated,
    delta_mean = mean_treated - mean_untreated,
    p_value = safe_wilcox(cur_data()),
    .groups="drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method="BH")) %>%
  arrange(suppressWarnings(as.numeric(cluster)), cluster)

# ------------------------------------------------------------------------------
# SAVE EXCEL
# ------------------------------------------------------------------------------
dir.create(file.path(outdir, "cell_level"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "sample_level"), showWarnings = FALSE, recursive = TRUE)

wb1 <- createWorkbook()
addWorksheet(wb1, "cell_stats")
addWorksheet(wb1, "cell_table")
writeData(wb1, "cell_stats", cell_stats)
writeData(wb1, "cell_table", df_cell %>% arrange(cluster, condition, sample))
saveWorkbook(wb1, file.path(outdir, "cell_level", "ISG_CELL_level.xlsx"), overwrite = TRUE)

wb2 <- createWorkbook()
addWorksheet(wb2, "sample_stats")
addWorksheet(wb2, "sample_table")
writeData(wb2, "sample_stats", sample_stats)
writeData(wb2, "sample_table", df_sample %>% arrange(cluster, condition, sample))
saveWorkbook(wb2, file.path(outdir, "sample_level", "ISG_SAMPLE_level.xlsx"), overwrite = TRUE)

# ------------------------------------------------------------------------------
# SAVE VIOLINS (CELL + SAMPLE)
# ------------------------------------------------------------------------------
cell_plot_dir <- file.path(outdir, "cell_level", "overlay_violins")
samp_plot_dir <- file.path(outdir, "sample_level", "overlay_violins")
dir.create(cell_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(samp_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Use ONE canvas width (A width) for both panels for identical layout
wA <- max(8.0, 0.55 * length(clusters_A) + 3.5)
wB <- wA

# ------------------------------------------------------------------------------
# CELL A
# ------------------------------------------------------------------------------
pA_cell <- plot_overlay(
  df_cell, clusters_A,
  title_text = "ISG Module Score across clusters per condition (CELL level) - panel A",
  subtitle_text = "Control (untreated) = grey; OT-I Treated = blue; close overlap; median stubs",
  ylab_text = sprintf("%s (cells; clipped 5-95%%)", score_col),
  nudge_val = nudge_A,
  violin_width = violin_width_A,
  median_half_width = median_half_A
)
ggsave(file.path(cell_plot_dir, "CELL_overlay_A.png"), pA_cell, width = wA, height = 5.0, dpi = 300)
ggsave(file.path(cell_plot_dir, "CELL_overlay_A.pdf"), pA_cell, width = wA, height = 5.0)

# ------------------------------------------------------------------------------
# CELL B (scaled geometry)
# ------------------------------------------------------------------------------
pB_cell <- plot_overlay(
  df_cell, clusters_B,
  title_text = "ISG Module Score across clusters per condition (CELL level) - panel B",
  subtitle_text = "Control (untreated) = grey; OT-I Treated = blue; close overlap; median stubs",
  ylab_text = sprintf("%s (cells; clipped 5-95%%)", score_col),
  nudge_val = nudge_B,
  violin_width = violin_width_B,
  median_half_width = median_half_B
)
ggsave(file.path(cell_plot_dir, "CELL_overlay_B.png"), pB_cell, width = wB, height = 5.0, dpi = 300)
ggsave(file.path(cell_plot_dir, "CELL_overlay_B.pdf"), pB_cell, width = wB, height = 5.0)

# ------------------------------------------------------------------------------
# SAMPLE A
# ------------------------------------------------------------------------------
pA_samp <- plot_overlay(
  df_sample, clusters_A,
  title_text = "ISG Module Score across clusters per condition (SAMPLE level) - panel A",
  subtitle_text = "Each sample = hashtag; value = median per (samplexcluster); Control grey; Treated blue",
  ylab_text = sprintf("%s (median per sample; clipped 5-95%%)", score_col),
  nudge_val = nudge_A,
  violin_width = violin_width_A,
  median_half_width = median_half_A
)
ggsave(file.path(samp_plot_dir, "SAMPLE_overlay_A.png"), pA_samp, width = wA, height = 5.0, dpi = 300)
ggsave(file.path(samp_plot_dir, "SAMPLE_overlay_A.pdf"), pA_samp, width = wA, height = 5.0)

# ------------------------------------------------------------------------------
# SAMPLE B (scaled geometry)
# ------------------------------------------------------------------------------
pB_samp <- plot_overlay(
  df_sample, clusters_B,
  title_text = "ISG Module Score across clusters per condition (SAMPLE level) - panel B",
  subtitle_text = "Each sample = hashtag; value = median per (samplexcluster); Control grey; Treated blue",
  ylab_text = sprintf("%s (median per sample; clipped 5-95%%)", score_col),
  nudge_val = nudge_B,
  violin_width = violin_width_B,
  median_half_width = median_half_B
)
ggsave(file.path(samp_plot_dir, "SAMPLE_overlay_B.png"), pB_samp, width = wB, height = 5.0, dpi = 300)
ggsave(file.path(samp_plot_dir, "SAMPLE_overlay_B.pdf"), pB_samp, width = wB, height = 5.0)

message("[OK] Done. Outputs saved under: ", normalizePath(outdir))
message(sprintf("Panel B geometry scale applied: %.4f (nB/nA = %d/%d)", scale_B, length(clusters_B), length(clusters_A)))
message("Violin fix applied: geom_violin(scale = 'width') in both conditions.")
message("Thin style applied: violin_width_A=0.70; nudge_A=0.04; median_half_A=0.07 (B scaled accordingly).")

