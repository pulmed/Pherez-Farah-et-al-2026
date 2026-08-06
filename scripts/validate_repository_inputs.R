#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/validate_repository_inputs.R
# Original file: repository validation utility
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Validate that the repository has the canonical Seurat object and required fields for downstream analyses.
# Inputs: Canonical analysis-ready Seurat object, default data/seurat.rds.
# Outputs: Console validation report with assay, layer, reduction, metadata, and count summaries.
# Assay/layer input: RNA counts/data, SCT data, ADT counts/data, optional HTO counts/data.
# Dependencies: Seurat, SeuratObject.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - This script performs validation only; it does not write analysis outputs.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("scripts/utils/seurat_io.R")

seurat_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
cat("Validating Seurat object:", seurat_path, "\n")

obj <- load_seurat_object(seurat_path)
summary <- object_summary(obj)

cat("Cells:", summary$cells, "\n")
cat("Features:", summary$features, "\n")
cat("Assays:", paste(summary$assays, collapse = ", "), "\n")
cat("Reductions:", paste(summary$reductions, collapse = ", "), "\n")

required_metadata <- c("final_clusters", "condition", "hash.ID")
required_assays <- c("RNA", "SCT", "ADT")
required_reductions <- c("pca", "umap")
required_layers <- list(
  RNA = c("counts", "data"),
  SCT = c("data"),
  ADT = c("counts", "data")
)

require_metadata(obj, required_metadata)
require_assays(obj, required_assays)
require_reductions(obj, required_reductions)
require_assay_layers(obj, required_layers)

cat("Required metadata: OK\n")
cat("Required assays: OK\n")
cat("Required reductions: OK\n")
cat("Required assay layers: OK\n")


cat("\nCluster counts:\n")
print(table(obj$final_clusters, useNA = "ifany"))

cat("\nCondition counts:\n")
print(table(obj$condition, useNA = "ifany"))

cat("\nSample counts:\n")
print(table(obj$hash.ID, useNA = "ifany"))

if ("HTO" %in% Assays(obj)) {
  cat("\nOptional HTO assay: present\n")
} else {
  cat("\nOptional HTO assay: absent\n")
}

cat("\nRepository input validation: OK\n")
