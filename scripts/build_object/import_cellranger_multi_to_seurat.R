#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/build_object/import_cellranger_multi_to_seurat.R
# Original file: process_seurat.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Import Cell Ranger multi outputs and merge samples into a multimodal Seurat object.
# Inputs: Cell Ranger multi output directories containing GEX, ADT/HTO, and optional V(D)J-T data.
# Outputs: Merged multimodal Seurat object.
# Assay/layer input: RNA counts, optional ADT counts, optional HTO counts, optional TCR metadata.
# Dependencies: Seurat, dplyr, scRepertoire, Matrix.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Override CELLRANGER_MULTI_DIR and OUTPUT_DIR as needed.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scRepertoire)
  library(Matrix)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------

# Root directory containing Cell Ranger multi outputs
# Expected structure:
# base_dir/
#   sample1_multi/
#   sample2_multi/
base_dir <- Sys.getenv("CELLRANGER_MULTI_DIR", unset = file.path("output", "cellranger_multi"))

# Output file
output_file <- file.path(
  Sys.getenv("OUTPUT_DIR", unset = "output"),
  "import_cellranger_multi_to_seurat",
  "merged_cellranger_multi_seurat.rds"
)

# Feature-detection settings
antibody_assay_name <- "ADT"
hto_assay_name <- "HTO"
hto_pattern <- "^Hashtag-"

# Whether to include TCR information when available
include_tcr <- TRUE

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

detect_samples <- function(base_dir) {
  sample_ids <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  sort(sample_ids[grepl("_multi$", sample_ids)])
}

get_matrix_path <- function(base_dir, sample_id) {
  file.path(
    base_dir,
    sample_id,
    "outs",
    "per_sample_outs",
    sample_id,
    "count",
    "sample_filtered_feature_bc_matrix"
  )
}

get_tcr_path <- function(base_dir, sample_id) {
  file.path(
    base_dir,
    sample_id,
    "outs",
    "per_sample_outs",
    sample_id,
    "vdj_t",
    "filtered_contig_annotations.csv"
  )
}

split_adt_and_hto <- function(antibody_capture, hto_pattern = "^Hashtag-") {
  hto_features <- grep(hto_pattern, rownames(antibody_capture), value = TRUE)
  adt_features <- setdiff(rownames(antibody_capture), hto_features)

  list(adt = adt_features, hto = hto_features)
}

add_optional_modalities <- function(
  seu,
  raw_data,
  hto_pattern,
  antibody_assay_name,
  hto_assay_name
) {
  if (!"Antibody Capture" %in% names(raw_data)) {
    message("No Antibody Capture assay found.")
    return(seu)
  }

  antibody_capture <- raw_data[["Antibody Capture"]]
  feature_split <- split_adt_and_hto(antibody_capture, hto_pattern)

  if (length(feature_split$adt) > 0) {
    seu[[antibody_assay_name]] <- CreateAssayObject(
      counts = antibody_capture[feature_split$adt, , drop = FALSE]
    )
  }

  if (length(feature_split$hto) > 0) {
    seu[[hto_assay_name]] <- CreateAssayObject(
      counts = antibody_capture[feature_split$hto, , drop = FALSE]
    )
  }

  seu
}

add_optional_tcr <- function(seu, tcr_path, sample_id) {
  if (!file.exists(tcr_path)) {
    message("No TCR file found for sample: ", sample_id)
    return(seu)
  }

  tcr <- read.csv(tcr_path, stringsAsFactors = FALSE)

  if (!"barcode" %in% colnames(tcr)) {
    warning("Missing 'barcode' column in TCR file for sample: ", sample_id)
    return(seu)
  }

  tcr_combined <- combineTCR(
    input.data = list(tcr),
    samples = sample_id
  )

  tcr_df <- tcr_combined[[1]]

  if (!all(c("barcode", "CTstrict") %in% colnames(tcr_df))) {
    warning("Invalid TCR structure for sample: ", sample_id)
    return(seu)
  }

  tcr_df$orig.ident <- tcr_df$barcode

  combineExpression(
    input.data = list(tcr_df),
    sc.data = seu,
    cloneCall = "CTstrict",
    group.by = "orig.ident"
  )
}

# ------------------------------------------------------------------------------
# 3. VALIDATION
# ------------------------------------------------------------------------------

if (!dir.exists(base_dir)) {
  stop("base_dir does not exist: ", base_dir)
}

sample_ids <- detect_samples(base_dir)

if (length(sample_ids) == 0) {
  stop("No '_multi' directories found in base_dir.")
}

message("Detected samples: ", paste(sample_ids, collapse = ", "))

# ------------------------------------------------------------------------------
# 4. IMPORT SAMPLES
# ------------------------------------------------------------------------------

seurat_list <- list()

for (sample_id in sample_ids) {
  message("Processing: ", sample_id)

  matrix_path <- get_matrix_path(base_dir, sample_id)

  if (!dir.exists(matrix_path)) {
    warning("Skipping (missing matrix): ", sample_id)
    next
  }

  raw_data <- Read10X(matrix_path)

  if (!"Gene Expression" %in% names(raw_data)) {
    warning("Skipping (no GEX): ", sample_id)
    next
  }

  seu <- CreateSeuratObject(
    counts = raw_data[["Gene Expression"]],
    project = sample_id
  )

  seu$sample_id <- sample_id

  seu <- add_optional_modalities(
    seu,
    raw_data,
    hto_pattern,
    antibody_assay_name,
    hto_assay_name
  )

  seu <- RenameCells(seu, add.cell.id = sample_id)

  if (include_tcr) {
    tcr_path <- get_tcr_path(base_dir, sample_id)
    seu <- add_optional_tcr(seu, tcr_path, sample_id)
  }

  seurat_list[[sample_id]] <- seu
}

# ------------------------------------------------------------------------------
# 5. MERGE AND SAVE
# ------------------------------------------------------------------------------

if (length(seurat_list) == 0) {
  stop("No Seurat objects created.")
}

combined_seurat <- if (length(seurat_list) == 1) {
  seurat_list[[1]]
} else {
  merge(seurat_list[[1]], y = seurat_list[-1])
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

saveRDS(combined_seurat, output_file)

message("Saved: ", output_file)
