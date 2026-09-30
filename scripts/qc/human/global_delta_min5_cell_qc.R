#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/qc/human/global_delta_min5_cell_qc.R
# Original file: 20260825 Global delta min 5 cells.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Build a human T0/T30 cluster-fraction delta QC table requiring minimum cells at both timepoints.
# Inputs: Human Seurat object with patient, timepoint, response, and cluster metadata.
# Outputs: Cell-count QC workbook, full/pass-fail delta tables, valid-patient summary, and waterfall plot.
# Assay/layer input: Metadata only.
# Dependencies: Seurat, dplyr, tidyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Result files are written under output/ and are not tracked by git.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
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
outdir <- make_output_dir(file.path(output_root, "human", "global_delta_min5_cell_qc"))

patient_col <- Sys.getenv("PATIENT_COL", unset = "patient_simple")
timepoint_col <- Sys.getenv("TIMEPOINT_COL", unset = "timepoint_simple")
cluster_col <- Sys.getenv("CLUSTER_COL", unset = "seurat_clusters")
response_col <- Sys.getenv("RESPONSE_COL", unset = "response_simple")

paired_patients <- parse_csv_env("PAIRED_PATIENTS", c("P1", "P5", "P6", "P10", "P7", "P8", "P13"))
timepoints_use <- parse_csv_env("TIMEPOINTS_USE", c("T0", "T30"))
excluded_clusters <- parse_csv_env("EXCLUDED_CLUSTERS", "6")
min_cells <- as.integer(Sys.getenv("MIN_CELLS", unset = "5"))

patient_colors <- c(
  P1 = "gray65",
  P10 = "gray45",
  P5 = "black",
  P6 = "gray20",
  P7 = "#9ECAE1",
  P8 = "#4292C6",
  P13 = "#08519C"
)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(rds_path)) {
  message("SMOKE_TEST=1: human Seurat object not found; skipping human global delta QC.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# LOAD METADATA
# ------------------------------------------------------------------------------
obj <- load_seurat_object(rds_path)
require_metadata(obj, c(patient_col, timepoint_col, cluster_col, response_col))

md <- obj@meta.data %>%
  as.data.frame() %>%
  transmute(
    patient = as.character(.data[[patient_col]]),
    timepoint = as.character(.data[[timepoint_col]]),
    cluster = as.character(.data[[cluster_col]]),
    response = as.character(.data[[response_col]])
  ) %>%
  filter(
    patient %in% paired_patients,
    timepoint %in% timepoints_use,
    !cluster %in% excluded_clusters
  )

patient_metadata <- md %>%
  distinct(patient, response) %>%
  arrange(response, patient)

clusters <- sort(unique(md$cluster))

# ------------------------------------------------------------------------------
# CELL COUNTS AND FRACTIONS
# ------------------------------------------------------------------------------
cluster_fractions <- md %>%
  count(patient, timepoint, cluster, name = "n_cells") %>%
  group_by(patient, timepoint) %>%
  mutate(total_cells = sum(n_cells), fraction = n_cells / total_cells) %>%
  ungroup()

complete_grid <- tidyr::expand_grid(
  patient = paired_patients,
  timepoint = timepoints_use,
  cluster = clusters
)

cluster_fractions <- complete_grid %>%
  left_join(cluster_fractions, by = c("patient", "timepoint", "cluster")) %>%
  group_by(patient, timepoint) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    total_cells = replace_na(total_cells, sum(n_cells)),
    fraction = ifelse(total_cells > 0, n_cells / total_cells, 0)
  ) %>%
  ungroup() %>%
  left_join(patient_metadata, by = "patient")

cell_count_wide <- cluster_fractions %>%
  mutate(patient_time = paste(patient, timepoint, sep = "_")) %>%
  select(cluster, patient_time, n_cells) %>%
  pivot_wider(names_from = patient_time, values_from = n_cells, values_fill = 0) %>%
  arrange(as.numeric(cluster))

delta_table <- cluster_fractions %>%
  select(patient, response, cluster, timepoint, n_cells, total_cells, fraction) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(n_cells, total_cells, fraction),
    values_fill = list(n_cells = 0, total_cells = 0, fraction = 0)
  ) %>%
  mutate(
    delta_T30_minus_T0 = fraction_T30 - fraction_T0,
    passes_min_cell_filter = n_cells_T0 >= min_cells & n_cells_T30 >= min_cells
  ) %>%
  arrange(as.numeric(cluster), response, patient)

excluded_under5 <- delta_table %>%
  filter(!passes_min_cell_filter)

valid_delta <- delta_table %>%
  filter(passes_min_cell_filter)

valid_patients <- valid_delta %>%
  count(cluster, response, name = "n_valid_patients") %>%
  arrange(as.numeric(cluster), response)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: global delta QC rows: ", nrow(delta_table))
  message("SMOKE_TEST=1: valid min-cell delta rows: ", nrow(valid_delta))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# EXPORT TABLES
# ------------------------------------------------------------------------------
workbook <- createWorkbook()
addWorksheet(workbook, "Cell counts")
writeData(workbook, "Cell counts", cell_count_wide)
addWorksheet(workbook, "Delta QC")
writeData(workbook, "Delta QC", delta_table)
addWorksheet(workbook, "Excluded under min")
writeData(workbook, "Excluded under min", excluded_under5)
addWorksheet(workbook, "Valid patients")
writeData(workbook, "Valid patients", valid_patients)
saveWorkbook(workbook, file.path(outdir, "GLOBAL_DELTA_MIN5_CELL_COUNT_QC.xlsx"), overwrite = TRUE)

write.csv(cell_count_wide, file.path(outdir, "CELL_COUNTS_PER_PATIENT_T0_T30.csv"), row.names = FALSE)
write.csv(delta_table, file.path(outdir, "GLOBAL_delta_MIN5_full_QC.csv"), row.names = FALSE)
write.csv(excluded_under5, file.path(outdir, "GLOBAL_delta_EXCLUDED_under5.csv"), row.names = FALSE)
write.csv(valid_patients, file.path(outdir, "VALID_PATIENTS_PER_CLUSTER_RESPONSE.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# WATERFALL PLOT
# ------------------------------------------------------------------------------
plot_df <- valid_delta %>%
  group_by(cluster) %>%
  arrange(delta_T30_minus_T0, .by_group = TRUE) %>%
  mutate(patient_order = factor(patient, levels = unique(patient))) %>%
  ungroup() %>%
  mutate(
    cluster = factor(cluster, levels = sort(unique(cluster))),
    response = factor(response, levels = c("Non_responder", "Responder"))
  )

p <- ggplot(plot_df, aes(x = patient_order, y = delta_T30_minus_T0, fill = response, color = patient)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey45") +
  geom_col(width = 0.72, alpha = 0.78, color = NA) +
  geom_point(size = 1.8) +
  facet_wrap(~ cluster, scales = "free_x") +
  scale_color_manual(values = patient_colors[names(patient_colors) %in% unique(plot_df$patient)], na.translate = FALSE) +
  labs(
    title = "Global delta composition with minimum cell-count QC",
    subtitle = paste0("Patient-cluster deltas require at least ", min_cells, " cells at both T0 and T30"),
    x = "Patient",
    y = "Delta fraction (T30 - T0)",
    fill = "Response",
    color = "Patient"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(outdir, "GLOBAL_delta_waterfall_MIN5_CELLS.pdf"), p, width = 12, height = 8)
ggsave(file.path(outdir, "GLOBAL_delta_waterfall_MIN5_CELLS.png"), p, width = 12, height = 8, dpi = 300)

message("Done. Outputs written to: ", outdir)
