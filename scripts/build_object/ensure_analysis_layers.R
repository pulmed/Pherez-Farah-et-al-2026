#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/build_object/ensure_analysis_layers.R
# Original file: repository utility helper
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Ensure the canonical Seurat object contains normalized assay layers expected by downstream analyses.
# Inputs: Canonical analysis-ready Seurat object, default data/seurat.rds.
# Outputs: Updated Seurat object with populated RNA data and, when present, normalized ADT/HTO data layers.
# Assay/layer input: RNA counts, optional ADT counts, optional HTO counts.
# Dependencies: Seurat, SeuratObject.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - By default this script updates SEURAT_RDS in place using a temporary file.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("scripts/utils/seurat_io.R")

input_file <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_file <- Sys.getenv("OUTPUT_RDS", unset = input_file)

data_layer_empty <- function(object, assay) {
  !assay_layer_exists(object, assay, "data")
}

if (!file.exists(input_file)) {
  stop("Input Seurat object does not exist: ", input_file, call. = FALSE)
}

obj <- readRDS(input_file)

if ("RNA" %in% Assays(obj) && data_layer_empty(obj, "RNA")) {
  message("Populating RNA:data with LogNormalize.")
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
} else {
  message("RNA:data already present; leaving unchanged.")
}

if ("ADT" %in% Assays(obj) && data_layer_empty(obj, "ADT")) {
  message("Populating ADT:data with CLR normalization.")
  obj <- NormalizeData(
    obj,
    assay = "ADT",
    normalization.method = "CLR",
    margin = 2,
    verbose = FALSE
  )
} else if ("ADT" %in% Assays(obj)) {
  message("ADT:data already present; leaving unchanged.")
}

if ("HTO" %in% Assays(obj) && data_layer_empty(obj, "HTO")) {
  message("Populating HTO:data with CLR normalization.")
  obj <- NormalizeData(
    obj,
    assay = "HTO",
    normalization.method = "CLR",
    margin = 2,
    verbose = FALSE
  )
} else if ("HTO" %in% Assays(obj)) {
  message("HTO:data already present; leaving unchanged.")
}

make_output_dir(dirname(output_file))
tmp_file <- paste0(output_file, ".tmp")
saveRDS(obj, tmp_file)
file.rename(tmp_file, output_file)

message("Saved updated Seurat object: ", output_file)
