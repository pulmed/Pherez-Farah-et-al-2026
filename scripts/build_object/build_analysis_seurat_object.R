#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/build_object/build_analysis_seurat_object.R
# Original file: normalization_seurat.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Build the canonical analysis-ready Seurat object with raw counts, normalized layers, SCT, reductions, and required metadata.
# Inputs: Merged multimodal Seurat object.
# Outputs: Canonical analysis-ready Seurat object, default data/seurat.rds, plus object-building QC plots.
# Assay/layer input: RNA counts/data, SCT data, ADT counts/data, optional HTO counts/data.
# Dependencies: Seurat, ggplot2, scRepertoire.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - HTO demultiplexing runs only when an HTO assay is present.
# - Batch correction is intentionally not applied because treatment may be confounded with sample.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scRepertoire)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------

# Input merged Seurat object
input_file <- Sys.getenv("MERGED_SEURAT_RDS", unset = file.path("data", "merged_cellranger_multi_seurat.rds"))

# Output normalized Seurat object
output_file <- Sys.getenv("SEURAT_RDS", unset = file.path("data", "seurat.rds"))

# Output directory for plots
plot_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path("output", "build_analysis_seurat_object"))

# Assay names
rna_assay <- "RNA"
adt_assay <- "ADT"
hto_assay <- "HTO"

# HTO demultiplexing
hto_positive_quantile <- 0.99
keep_only_hto_singlets <- TRUE

# Condition assignment from HTO labels
treated_hashtags <- paste0("Hashtag-0", 1:4, "-TotalSeqC")
treated_label <- "treated"
untreated_label <- "untreated"

# QC thresholds
min_features_rna <- 200
max_percent_mt <- 10
mitochondrial_gene_pattern <- "^mt-"

# RNA log-normalization
run_rna_log_normalization <- TRUE
rna_normalization_method <- "LogNormalize"
rna_scale_factor <- 10000

# Dimensionality reduction / clustering
dims_use <- 1:20
cluster_resolution <- 0.8

# Biotin-positive ADT threshold
biotin_feature_pattern <- "Biotin"
biotin_positive_threshold <- 1

# Feature plots
rna_features_to_plot <- c("Cd4", "Ifng")
adt_features_to_plot <- c("Biotin-TotalSeqC")

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height
  )
}

has_assay <- function(object, assay_name) {
  assay_name %in% Assays(object)
}

print_object_summary <- function(object) {
  message("Assays:")
  print(Assays(object))

  message("Metadata columns:")
  print(colnames(object@meta.data))

  message("Number of cells: ", ncol(object))
  message("Default assay: ", DefaultAssay(object))
}

assign_condition <- function(object, treated_hashtags, treated_label, untreated_label) {
  object$sample_id <- object$HTO_classification
  object$condition <- ifelse(
    object$sample_id %in% treated_hashtags,
    treated_label,
    untreated_label
  )

  object
}

count_biotin_positive_cells <- function(
  object,
  adt_assay,
  biotin_feature_pattern,
  biotin_positive_threshold
) {
  if (!has_assay(object, adt_assay)) {
    message("ADT assay not found; skipping biotin-positive cell count.")
    return(invisible(NULL))
  }

  biotin_features <- grep(
    biotin_feature_pattern,
    rownames(object[[adt_assay]]),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(biotin_features) == 0) {
    message("No biotin-labelled features found in ADT assay.")
    return(invisible(NULL))
  }

  adt_data <- GetAssayData(object, assay = adt_assay, slot = "data")

  for (feature in biotin_features) {
    n_positive <- sum(adt_data[feature, ] > biotin_positive_threshold)

    message(
      "Biotin-positive cells for ",
      feature,
      " (ADT > ",
      biotin_positive_threshold,
      "): ",
      n_positive
    )
  }

  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 3. LOAD DATA AND VALIDATE INPUTS
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

combined <- readRDS(input_file)

if (!has_assay(combined, rna_assay)) {
  stop("RNA assay not found: ", rna_assay)
}

DefaultAssay(combined) <- rna_assay

print_object_summary(combined)

# ------------------------------------------------------------------------------
# 4. HTO DEMULTIPLEXING
# ------------------------------------------------------------------------------

if (has_assay(combined, hto_assay)) {
  message("Running HTO normalization and demultiplexing.")

  combined <- NormalizeData(
    combined,
    assay = hto_assay,
    normalization.method = "CLR"
  )

  combined <- HTODemux(
    combined,
    assay = hto_assay,
    positive.quantile = hto_positive_quantile
  )

  message("HTO demultiplexing result:")
  print(table(combined$HTO_classification.global))

  if (keep_only_hto_singlets) {
    combined <- subset(
      combined,
      subset = HTO_classification.global == "Singlet"
    )
    message("Retained HTO singlets: ", ncol(combined))
  }

  combined <- assign_condition(
    object = combined,
    treated_hashtags = treated_hashtags,
    treated_label = treated_label,
    untreated_label = untreated_label
  )
} else {
  message("HTO assay not found; skipping demultiplexing.")
}

# ------------------------------------------------------------------------------
# 5. QUALITY CONTROL AND NORMALIZATION
# ------------------------------------------------------------------------------

DefaultAssay(combined) <- rna_assay

if (!"percent.mt" %in% colnames(combined@meta.data)) {
  combined[["percent.mt"]] <- PercentageFeatureSet(
    combined,
    pattern = mitochondrial_gene_pattern
  )
}

combined <- subset(
  combined,
  subset = nFeature_RNA > min_features_rna & percent.mt < max_percent_mt
)

message("Cells after QC filtering: ", ncol(combined))

# Store standard log-normalized RNA values for downstream plotting/scripts.
# This fills the RNA assay data slot and avoids re-running NormalizeData later.
if (run_rna_log_normalization) {
  DefaultAssay(combined) <- rna_assay

  combined <- NormalizeData(
    combined,
    assay = rna_assay,
    normalization.method = rna_normalization_method,
    scale.factor = rna_scale_factor,
    verbose = FALSE
  )
}

# Normalize ADT data if available.
if (has_assay(combined, adt_assay)) {
  combined <- NormalizeData(
    combined,
    assay = adt_assay,
    normalization.method = "CLR",
    margin = 2,
    verbose = FALSE
  )
}

# SCTransform is used for dimensionality reduction, clustering, and UMAP.
combined <- SCTransform(
  combined,
  assay = rna_assay,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

DefaultAssay(combined) <- "SCT"

# ------------------------------------------------------------------------------
# 6. BIOTIN QC
# ------------------------------------------------------------------------------

count_biotin_positive_cells(
  object = combined,
  adt_assay = adt_assay,
  biotin_feature_pattern = biotin_feature_pattern,
  biotin_positive_threshold = biotin_positive_threshold
)

# ------------------------------------------------------------------------------
# 7. DIMENSIONALITY REDUCTION, CLUSTERING, AND PLOTS
# ------------------------------------------------------------------------------

combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, dims = dims_use, verbose = FALSE)
combined <- FindNeighbors(combined, dims = dims_use, verbose = FALSE)
combined <- FindClusters(
  combined,
  resolution = cluster_resolution,
  verbose = FALSE
)

if ("condition" %in% colnames(combined@meta.data)) {
  p_condition <- DimPlot(combined, group.by = "condition")
  save_plot(p_condition, "umap_by_condition.png")

  p_clusters_split <- DimPlot(
    combined,
    split.by = "condition",
    group.by = "seurat_clusters"
  )
  save_plot(p_clusters_split, "umap_clusters_split_by_condition.png", width = 10)
}

DefaultAssay(combined) <- rna_assay

available_rna_features <- intersect(
  rna_features_to_plot,
  rownames(combined[[rna_assay]])
)

if (length(available_rna_features) > 0) {
  p_rna_features <- FeaturePlot(
    combined,
    features = available_rna_features,
    split.by = if ("condition" %in% colnames(combined@meta.data)) "condition" else NULL,
    max.cutoff = "q95"
  )

  save_plot(p_rna_features, "rna_feature_plots.png", width = 10)
} else {
  message("None of the configured RNA features were found.")
}

if (has_assay(combined, adt_assay)) {
  DefaultAssay(combined) <- adt_assay

  available_adt_features <- intersect(
    adt_features_to_plot,
    rownames(combined[[adt_assay]])
  )

  if (length(available_adt_features) > 0) {
    p_adt_features <- FeaturePlot(
      combined,
      features = available_adt_features,
      split.by = if ("condition" %in% colnames(combined@meta.data)) "condition" else NULL,
      min.cutoff = "q05",
      max.cutoff = "q95"
    )

    save_plot(p_adt_features, "adt_feature_plots.png", width = 10)
  } else {
    message("None of the configured ADT features were found.")
  }
}

DefaultAssay(combined) <- "SCT"

# ------------------------------------------------------------------------------
# 8. SAVE OUTPUT
# ------------------------------------------------------------------------------

saveRDS(combined, file = output_file)

message("Saved normalized Seurat object to: ", output_file)
