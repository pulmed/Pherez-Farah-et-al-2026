#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/cross_species/mouse_to_human/score_cluster2_isg_responder_status.R
# Original file: 20260825 cluster 2 ISG.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Compare human cluster-2 ISG scores between responders and non-responders.
# Inputs: Human Seurat object with RNA data and patient/timepoint/response/cluster metadata.
# Outputs: Patient-level ISG workbook and responder-vs-non-responder plots.
# Assay/layer input: RNA data.
# Dependencies: Seurat, SeuratObject, dplyr, tidyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - The all-timepoint cluster-2 comparison is the primary manuscript comparison.
# - Uses the final 47-gene high-confidence human ISG ortholog set.
# - This script lives under mouse_to_human because it interrogates an IFN/ISG program
#   in the human data for comparison with the mouse signature results.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

source("scripts/utils/seurat_io.R")
source("scripts/utils/signature_helpers.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv(
  "HUMAN_ANALYSIS_RDS",
  unset = Sys.getenv("HUMAN_RPCA_RDS", unset = "data/human_seurat_rpca.rds")
)
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "human", "cluster2_isg_responder_status"))

assay_use <- Sys.getenv("ASSAY", unset = "RNA")
patient_col <- Sys.getenv("PATIENT_COL", unset = "patient_simple")
timepoint_col <- Sys.getenv("TIMEPOINT_COL", unset = "timepoint_simple")
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "seurat_clusters")
response_col <- Sys.getenv("RESPONSE_COL", unset = "response_simple")
cluster_use <- Sys.getenv("CLUSTER_USE", unset = "2")
patients <- parse_csv_env("PATIENTS", c("P1", "P5", "P6", "P10", "P7", "P8", "P13"))
timepoints_use <- parse_csv_env("TIMEPOINTS_USE", c("T0", "T30"))
response_levels <- parse_csv_env("RESPONSE_LEVELS", c("Non_responder", "Responder"))
score_seed <- as.integer(Sys.getenv("MODULE_SCORE_SEED", unset = "1234"))

patient_colors <- c(
  P1 = "gray65",
  P10 = "gray45",
  P5 = "black",
  P6 = "gray20",
  P7 = "#9ECAE1",
  P8 = "#4292C6",
  P13 = "#08519C"
)

isg_genes <- c(
  "ARID5B", "ATP11B", "ATP6V0C", "BACH1", "BATF2", "CD300LF",
  "CEBPB", "CREB5", "CTSS", "CXCL10", "CYRIB", "ENTPD1",
  "FCGR1A", "FCGR3B", "FCGR3A", "FGL2", "GBP2", "GNA13",
  "GPR141", "PYHIN1", "MNDA", "IFI27L2", "IFIT2", "IRGM",
  "IL10RA", "IRAK2", "IRF7", "LILRB4", "OAS3", "PARP14",
  "PDE7B", "PIK3AP1", "PNP", "PTPRC", "RNF213", "RTP4",
  "SAMD9L", "SAMHD1", "SDCBP", "SLAMF8", "SLFN13", "STAT2",
  "TGFBI", "TGM2", "TNFAIP2", "XAF1", "ZBP1"
)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(rds_path)) {
  message("SMOKE_TEST=1: human Seurat object not found; skipping human cluster 2 ISG scoring.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# LOAD AND SCORE
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_assays(obj, assay_use)
require_metadata(obj, c(patient_col, timepoint_col, cluster_col, response_col))

DefaultAssay(obj) <- assay_use
obj <- coerce_metadata_character(obj, c(patient_col, timepoint_col, cluster_col, response_col))
obj <- ensure_rna_data_layer(obj, assay = assay_use)

features <- rownames(obj[[assay_use]])
present_genes <- present_signature_genes(isg_genes, features, "Cluster2_ISG")
gene_check <- gene_presence_table(isg_genes, features, "Cluster2_ISG")

score <- add_module_score(obj, present_genes, "ISG", assay = assay_use, seed = score_seed)
obj <- score$object

df <- obj@meta.data %>%
  as.data.frame() %>%
  transmute(
    patient = as.character(.data[[patient_col]]),
    response = as.character(.data[[response_col]]),
    timepoint = as.character(.data[[timepoint_col]]),
    cluster = as.character(.data[[cluster_col]]),
    ISG = .data[[score$column]]
  ) %>%
  filter(
    patient %in% patients,
    cluster == cluster_use,
    timepoint %in% timepoints_use,
    response %in% response_levels
  )

if (nrow(df) == 0) {
  stop("No cells remained after cluster/patient/timepoint filtering.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# SUMMARIES AND STATS
# ------------------------------------------------------------------------------
make_stats <- function(patient_df, analysis_name) {
  responder <- patient_df$mean_ISG[patient_df$response == "Responder"]
  non_responder <- patient_df$mean_ISG[patient_df$response == "Non_responder"]

  data.frame(
    analysis = analysis_name,
    n_responders = length(responder),
    n_non_responders = length(non_responder),
    mean_responder = mean(responder, na.rm = TRUE),
    mean_non_responder = mean(non_responder, na.rm = TRUE),
    median_responder = median(responder, na.rm = TRUE),
    median_non_responder = median(non_responder, na.rm = TRUE),
    delta_mean_responder_minus_non_responder = mean(responder, na.rm = TRUE) - mean(non_responder, na.rm = TRUE),
    test = "Wilcoxon rank-sum",
    p_value = safe_wilcox_p(responder, non_responder)
  )
}

all_timepoints_patient <- df %>%
  group_by(patient, response) %>%
  summarise(
    mean_ISG = mean(ISG, na.rm = TRUE),
    median_ISG = median(ISG, na.rm = TRUE),
    n_cells = n(),
    n_T0_cells = sum(timepoint == "T0"),
    n_T30_cells = sum(timepoint == "T30"),
    .groups = "drop"
  )

t30_patient <- df %>%
  filter(timepoint == "T30") %>%
  group_by(patient, response) %>%
  summarise(
    mean_ISG = mean(ISG, na.rm = TRUE),
    median_ISG = median(ISG, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

all_timepoints_stats <- make_stats(all_timepoints_patient, "Cluster 2 ISG: responder vs non-responder, all timepoints")
t30_stats <- make_stats(t30_patient, "Cluster 2 ISG: responder vs non-responder, T30 only")

cell_count_qc <- df %>%
  count(patient, response, timepoint, name = "n_cells") %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = n_cells, values_fill = 0) %>%
  arrange(response, patient)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: cluster 2 ISG genes present: ", length(present_genes), "/", length(isg_genes))
  message("SMOKE_TEST=1: filtered cluster 2 cells: ", nrow(df))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------------------
make_plot <- function(patient_df, title_text) {
  plot_df <- patient_df %>%
    mutate(response = factor(response, levels = response_levels))

  ggplot(plot_df, aes(x = response, y = mean_ISG, color = patient)) +
    geom_boxplot(aes(group = response), color = "grey35", fill = "white", outlier.shape = NA, width = 0.45) +
    geom_point(size = 3) +
    geom_text(aes(label = patient), vjust = -0.8, size = 3) +
    scale_color_manual(values = patient_colors[names(patient_colors) %in% unique(plot_df$patient)], na.translate = FALSE) +
    labs(title = title_text, x = NULL, y = "Mean ISG module score per patient", color = "Patient") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold"))
}

p_all <- make_plot(all_timepoints_patient, "Cluster 2 ISG: responder vs non-responder, all timepoints")
p_t30 <- make_plot(t30_patient, "Cluster 2 ISG: responder vs non-responder, T30 only")

ggsave(file.path(outdir, "CLUSTER2_ISG_R_vs_NR_ALL_TIMEPOINTS.pdf"), p_all, width = 5.5, height = 5)
ggsave(file.path(outdir, "CLUSTER2_ISG_R_vs_NR_ALL_TIMEPOINTS.png"), p_all, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(outdir, "CLUSTER2_ISG_R_vs_NR_T30_ONLY.pdf"), p_t30, width = 5.5, height = 5)
ggsave(file.path(outdir, "CLUSTER2_ISG_R_vs_NR_T30_ONLY.png"), p_t30, width = 5.5, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# EXPORTS
# ------------------------------------------------------------------------------
workbook <- createWorkbook()
addWorksheet(workbook, "All T patient values")
writeData(workbook, "All T patient values", all_timepoints_patient)
addWorksheet(workbook, "T30 patient values")
writeData(workbook, "T30 patient values", t30_patient)
addWorksheet(workbook, "Primary all T stats")
writeData(workbook, "Primary all T stats", all_timepoints_stats)
addWorksheet(workbook, "T30 stats")
writeData(workbook, "T30 stats", t30_stats)
addWorksheet(workbook, "Cell count QC")
writeData(workbook, "Cell count QC", cell_count_qc)
addWorksheet(workbook, "Gene presence")
writeData(workbook, "Gene presence", gene_check)
saveWorkbook(workbook, file.path(outdir, "CLUSTER2_ISG_COMPLETE_ANALYSIS.xlsx"), overwrite = TRUE)

message("Done. Outputs written to: ", outdir)
