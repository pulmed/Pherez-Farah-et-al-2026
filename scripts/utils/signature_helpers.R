#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/utils/signature_helpers.R
# Original file: repository utility helper
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Provide shared helpers for signature-scoring scripts.
# Inputs: Seurat objects, signature gene vectors, metadata columns, and output paths.
# Outputs: Gene-presence tables, module-score columns, and statistical summaries.
# Assay/layer input: Helper-dependent; callers choose the expression assay/layer.
# Dependencies: Seurat, SeuratObject.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - This file is intended to be sourced by scripts, not run directly.
# ------------------------------------------------------------------------------

parse_csv_env <- function(name, default) {
  value <- Sys.getenv(name, unset = paste(default, collapse = ","))
  value <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  value[nzchar(value)]
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

first_existing_metadata <- function(object, candidates) {
  hit <- candidates[candidates %in% colnames(object@meta.data)]
  if (length(hit) == 0) {
    stop("Missing required metadata column; tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }

  hit[[1]]
}

coerce_metadata_character <- function(object, columns) {
  for (column in columns) {
    object[[column]] <- as.character(object[[column]][, 1])
  }

  object
}

gene_presence_table <- function(genes, features, signature = NULL) {
  present <- genes %in% features
  out <- data.frame(
    gene = genes,
    status = ifelse(present, "Present", "Missing"),
    stringsAsFactors = FALSE
  )

  if (!is.null(signature)) {
    out <- data.frame(signature = signature, out, stringsAsFactors = FALSE)
  }

  out
}

present_signature_genes <- function(genes, features, signature) {
  present <- intersect(genes, features)
  if (length(present) == 0) {
    stop("No genes from signature found in expression matrix: ", signature, call. = FALSE)
  }

  present
}

drop_metadata_column <- function(object, column) {
  if (column %in% colnames(object@meta.data)) {
    object[[column]] <- NULL
  }

  object
}

add_module_score <- function(object, genes, name, assay = "RNA", seed = 1234) {
  score_col <- paste0(name, "1")
  object <- drop_metadata_column(object, score_col)

  set.seed(seed)
  object <- Seurat::AddModuleScore(
    object = object,
    features = list(genes),
    assay = assay,
    name = name,
    seed = seed,
    verbose = FALSE
  )

  list(object = object, column = score_col)
}

safe_wilcox_p <- function(x, y, paired = FALSE) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }

  if (paired && length(x) != length(y)) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(x, y, paired = paired, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

p_stars <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.0001 ~ "****",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}
