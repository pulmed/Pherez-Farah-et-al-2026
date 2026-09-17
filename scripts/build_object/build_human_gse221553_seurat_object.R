#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/build_object/build_human_gse221553_seurat_object.R
# Original file: GSE221553_Seurat.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Build a local human Seurat object from GSE221553 per-sample count tables.
# Inputs: GEO count tables named *-counts.tsv.gz under GSE221553_DIR, optionally downloaded from GEO.
# Outputs: Local human Seurat RDS and optional QC plots.
# Assay/layer input: RNA counts.
# Dependencies: Seurat, SeuratObject, Matrix, ggplot2.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - The generated human object is local derived data and is not tracked by git.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

source("scripts/utils/seurat_io.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
data_dir <- Sys.getenv("GSE221553_DIR", unset = "data/human/GSE221553")
output_rds <- Sys.getenv("HUMAN_SEURAT_RDS", unset = "data/human_seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "human", "build_gse221553"))

geo_accession <- Sys.getenv("GEO_ACCESSION", unset = "GSE221553")
download_geo <- identical(tolower(Sys.getenv("DOWNLOAD_GEO", unset = "false")), "true")
project_name <- Sys.getenv("PROJECT_NAME", unset = "GSE221553")
min_cells <- as.integer(Sys.getenv("MIN_CELLS", unset = "3"))
min_features_create <- as.integer(Sys.getenv("MIN_FEATURES_CREATE", unset = "200"))
min_features_qc <- as.integer(Sys.getenv("MIN_FEATURES_QC", unset = "200"))
max_features_qc <- as.integer(Sys.getenv("MAX_FEATURES_QC", unset = "7000"))
max_percent_mt <- as.numeric(Sys.getenv("MAX_PERCENT_MT", unset = "20"))
variable_features <- as.integer(Sys.getenv("VARIABLE_FEATURES", unset = "3000"))
dims_use <- seq_len(as.integer(Sys.getenv("DIMS", unset = "30")))
cluster_resolution <- as.numeric(Sys.getenv("CLUSTER_RESOLUTION", unset = "0.5"))
run_sct <- identical(tolower(Sys.getenv("RUN_SCT", unset = "true")), "true")
save_qc_plots <- identical(tolower(Sys.getenv("SAVE_QC_PLOTS", unset = "true")), "true")

count_files <- list.files(
  data_dir,
  pattern = "-counts[.]tsv[.]gz$",
  full.names = TRUE
)

geo_stub <- function(accession) {
  sub("[0-9]{1,3}$", "nnn", accession)
}

series_geo_url <- function(accession, section, filename = NULL) {
  url <- file.path(
    "https://ftp.ncbi.nlm.nih.gov/geo/series",
    geo_stub(accession),
    accession,
    section
  )

  if (!is.null(filename)) {
    url <- file.path(url, filename)
  }

  url
}

download_if_missing <- function(url, destfile) {
  if (file.exists(destfile)) {
    message("Using cached GEO file: ", destfile)
    return(invisible(destfile))
  }

  message("Downloading ", url)
  utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
  invisible(destfile)
}

download_geo_counts <- function(accession, destination_dir) {
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)

  soft_file <- paste0(accession, "_family.soft.gz")
  soft_path <- file.path(destination_dir, soft_file)
  soft_url <- series_geo_url(accession, "soft", soft_file)
  download_if_missing(soft_url, soft_path)

  soft_lines <- readLines(gzfile(soft_path), warn = FALSE)
  supp_lines <- grep("^!Sample_supplementary_file", soft_lines, value = TRUE)
  supp_urls <- sub("^!Sample_supplementary_file(_[0-9]+)? = ", "", supp_lines)
  supp_urls <- gsub("^ftp://ftp[.]ncbi[.]nlm[.]nih[.]gov", "https://ftp.ncbi.nlm.nih.gov", supp_urls)
  count_urls <- unique(supp_urls[grepl("-counts[.]tsv[.]gz$", supp_urls)])

  if (length(count_urls) == 0) {
    stop("No *-counts.tsv.gz supplementary files found in GEO SOFT file: ", soft_path, call. = FALSE)
  }

  for (url in count_urls) {
    destfile <- file.path(destination_dir, utils::URLdecode(basename(url)))
    download_if_missing(url, destfile)
  }

  count_urls
}

if (length(count_files) == 0 && download_geo) {
  download_geo_counts(geo_accession, data_dir)
  count_files <- list.files(
    data_dir,
    pattern = "-counts[.]tsv[.]gz$",
    full.names = TRUE
  )
}

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: GSE221553 count files found: ", length(count_files))
  quit(save = "no", status = 0)
}

if (length(count_files) == 0) {
  stop(
    "No GSE221553 count files found under ", data_dir,
    ". Expected files named *-counts.tsv.gz. Set GSE221553_DIR to override.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# LOAD PER-SAMPLE COUNT TABLES
# ------------------------------------------------------------------------------
extract_match <- function(x, pattern) {
  out <- sub(pattern, "\\1", x)
  ifelse(identical(out, x), NA_character_, out)
}

read_sample_object <- function(path) {
  sample_id <- sub("-counts[.]tsv[.]gz$", "", basename(path))
  message("Reading ", sample_id)

  mat <- read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)

  object <- CreateSeuratObject(
    counts = mat,
    project = project_name,
    min.cells = min_cells,
    min.features = min_features_create
  )

  object$sample_id <- sample_id
  object$gsm <- extract_match(sample_id, "^([^_]+)_.*$")
  object$patient <- extract_match(sample_id, ".*_(patient[0-9]+)-.*$")
  object$timepoint <- extract_match(sample_id, ".*-(T[0-9]+)-.*$")
  object
}

sample_objs <- lapply(count_files, read_sample_object)
names(sample_objs) <- sub("-counts[.]tsv[.]gz$", "", basename(count_files))

if (length(sample_objs) == 1) {
  obj <- sample_objs[[1]]
} else {
  obj <- merge(
    x = sample_objs[[1]],
    y = sample_objs[-1],
    add.cell.ids = names(sample_objs),
    project = project_name
  )
}

if (length(SeuratObject::Layers(obj[["RNA"]])) > 1) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

# ------------------------------------------------------------------------------
# QC AND NORMALIZATION
# ------------------------------------------------------------------------------
mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
if (length(mt_genes) > 0) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
} else {
  obj$percent.mt <- 0
  warning("No mitochondrial genes matched '^MT-'. percent.mt set to 0.")
}

if (save_qc_plots) {
  p_vln <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    layer = "counts"
  )
  ggsave(file.path(outdir, "GSE221553_QC_violin_prefilter.pdf"), p_vln, width = 10, height = 4.5)

  p_counts_features <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(file.path(outdir, "GSE221553_QC_counts_vs_features.pdf"), p_counts_features, width = 5.5, height = 5)

  p_counts_mt <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(file.path(outdir, "GSE221553_QC_counts_vs_percent_mt.pdf"), p_counts_mt, width = 5.5, height = 5)
}

obj <- subset(
  obj,
  subset = nFeature_RNA > min_features_qc &
    nFeature_RNA < max_features_qc &
    percent.mt < max_percent_mt
)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = variable_features, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

if (run_sct) {
  obj <- SCTransform(
    obj,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )

  obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, reduction = "pca", verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, reduction = "pca", verbose = FALSE)
}

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_rds)

summary <- data.frame(
  metric = c("cells", "features", "samples", "output_rds"),
  value = c(ncol(obj), nrow(obj), length(sample_objs), output_rds)
)
write.csv(summary, file.path(outdir, "GSE221553_build_summary.csv"), row.names = FALSE)

message("Done. Human object written to: ", output_rds)
