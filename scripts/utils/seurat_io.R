#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/utils/seurat_io.R
# Original file: repository utility helper
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Provide shared Seurat input/output, validation, and directory helpers.
# Inputs: Canonical analysis-ready Seurat object paths and required field names.
# Outputs: Loaded objects, validation errors, and created output directories.
# Assay/layer input: Helper-dependent; no assay is read unless requested by caller.
# Dependencies: Seurat.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - This file is intended to be sourced by scripts, not run directly.
# ------------------------------------------------------------------------------

load_seurat_object <- function(path = Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")) {
  if (!file.exists(path)) {
    stop("Seurat object not found: ", path, call. = FALSE)
  }

  readRDS(path)
}

make_output_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

require_metadata <- function(object, columns) {
  missing_columns <- setdiff(columns, colnames(object@meta.data))
  if (length(missing_columns) > 0) {
    stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

require_assays <- function(object, assays) {
  missing_assays <- setdiff(assays, Seurat::Assays(object))
  if (length(missing_assays) > 0) {
    stop("Missing required assays: ", paste(missing_assays, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

require_reductions <- function(object, reductions) {
  missing_reductions <- setdiff(reductions, Seurat::Reductions(object))
  if (length(missing_reductions) > 0) {
    stop("Missing required reductions: ", paste(missing_reductions, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

assay_layer_exists <- function(object, assay, layer) {
  if (!assay %in% Seurat::Assays(object)) return(FALSE)

  value <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )

  !is.null(value) && nrow(value) > 0 && ncol(value) > 0
}

require_assay_layers <- function(object, requirements) {
  missing_layers <- character()

  for (assay in names(requirements)) {
    for (layer in requirements[[assay]]) {
      if (!assay_layer_exists(object, assay, layer)) {
        missing_layers <- c(missing_layers, paste0(assay, ":", layer))
      }
    }
  }

  if (length(missing_layers) > 0) {
    stop("Missing required assay layers: ", paste(missing_layers, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

ensure_rna_data_layer <- function(object, assay = "RNA") {
  require_assays(object, assay)

  if (!assay_layer_exists(object, assay, "data")) {
    message("RNA data layer missing or empty; running NormalizeData().")
    object <- Seurat::NormalizeData(
      object,
      assay = assay,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }

  object
}

object_summary <- function(object) {
  list(
    cells = ncol(object),
    features = nrow(object),
    assays = Seurat::Assays(object),
    reductions = Seurat::Reductions(object),
    metadata_columns = colnames(object@meta.data)
  )
}
