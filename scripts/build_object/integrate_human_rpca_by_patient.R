#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/build_object/integrate_human_rpca_by_patient.R
# Original file: human RPCA processing provenance
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Perform patient-level RPCA integration of the local human GSE221553 Seurat object.
# Inputs: Human Seurat object with RNA assay and patient metadata.
# Outputs: Human Seurat object with integrated assay, PCA, UMAP, neighbors, and clusters.
# Assay/layer input: RNA counts/data.
# Dependencies: Seurat, SeuratObject.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Generated human objects are local derived data and are not tracked by git.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source("scripts/utils/seurat_io.R")
source("scripts/utils/signature_helpers.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
input_rds <- Sys.getenv("HUMAN_SEURAT_RDS", unset = "data/human_seurat.rds")
output_rds <- Sys.getenv("HUMAN_RPCA_RDS", unset = "data/human_seurat_rpca.rds")
patient_col_candidates <- parse_csv_env("PATIENT_COL_CANDIDATES", c("patient_simple", "patient", "patient_id"))
dims_use <- seq_len(as.integer(Sys.getenv("DIMS", unset = "30")))
cluster_resolution <- as.numeric(Sys.getenv("CLUSTER_RESOLUTION", unset = "0.5"))
nfeatures <- as.integer(Sys.getenv("VARIABLE_FEATURES", unset = "3000"))

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(input_rds)) {
  message("SMOKE_TEST=1: human Seurat object not found; skipping RPCA integration.")
  quit(save = "no", status = 0)
}

obj <- load_seurat_object(input_rds)
require_assays(obj, "RNA")
patient_col <- first_existing_metadata(obj, patient_col_candidates)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: human RPCA patient column selected: ", patient_col)
  message("SMOKE_TEST=1: cells: ", ncol(obj))
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# RPCA INTEGRATION
# ------------------------------------------------------------------------------
DefaultAssay(obj) <- "RNA"
obj <- coerce_metadata_character(obj, patient_col)
obj <- ensure_rna_data_layer(obj, assay = "RNA")

if (length(SeuratObject::Layers(obj[["RNA"]])) > 1) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

patient_objects <- SplitObject(obj, split.by = patient_col)
patient_objects <- lapply(patient_objects, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(object.list = patient_objects, nfeatures = nfeatures)
patient_objects <- lapply(patient_objects, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})

anchors <- FindIntegrationAnchors(
  object.list = patient_objects,
  anchor.features = features,
  reduction = "rpca",
  dims = dims_use,
  verbose = FALSE
)

obj_integrated <- IntegrateData(anchorset = anchors, dims = dims_use, verbose = FALSE)
DefaultAssay(obj_integrated) <- "integrated"
obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
obj_integrated <- RunPCA(obj_integrated, npcs = max(dims_use), verbose = FALSE)
obj_integrated <- FindNeighbors(obj_integrated, dims = dims_use, verbose = FALSE)
obj_integrated <- FindClusters(obj_integrated, resolution = cluster_resolution, verbose = FALSE)
obj_integrated <- RunUMAP(obj_integrated, dims = dims_use, verbose = FALSE)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_integrated, output_rds)

message("Done. RPCA-integrated human object written to: ", output_rds)
