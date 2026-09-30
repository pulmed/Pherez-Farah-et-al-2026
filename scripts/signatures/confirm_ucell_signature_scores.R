#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/confirm_ucell_signature_scores.R
# Original files: 20260928 * UCell.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Run UCell scoring as a confirmation analysis for selected final signatures.
# Inputs: Canonical mouse Seurat object and, optionally, a human Seurat object.
# Outputs: UCell gene-presence tables plus cell/sample-level score summaries.
# Assay/layer input: RNA counts when available; RNA data otherwise.
# Dependencies: Seurat, UCell, dplyr, tidyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Primary manuscript signature analyses use AddModuleScore; this script is a sensitivity/confirmation analysis.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

mouse_rds <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
human_rds <- Sys.getenv(
  "HUMAN_ANALYSIS_RDS",
  unset = Sys.getenv("HUMAN_RPCA_RDS", unset = "data/human_seurat_rpca.rds")
)
output_dir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "18_ucell_signature_confirmation")

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
condition_col <- "condition"
assay_name <- "RNA"

mouse_signatures <- list(
  ISG70 = c(
    "Irgm2", "Fcgr3", "Il10ra", "Pik3ap1", "Oas3", "Ifi47", "Entpd1", "Tnfaip2", "Rtp4", "Gna13",
    "Irgm1", "Fcgr4", "Stat2", "Sdcbp", "Gm4951", "Ifi213", "Samd9l", "Igtp", "Ifi207", "Xaf1",
    "Socs3", "Mndal", "Pde7b", "Irf7", "Gm12185", "Slamf8", "Atp11b", "Samhd1", "Fbxl5", "Tgtp2",
    "Slfn1", "Cd86", "Ifi27l2a", "Gm20663", "Zbp1", "Slfn8", "Creb5", "Irak2", "Ifi206", "Pnp",
    "Iigp1", "Ly6i", "Batf2", "Arid5b", "Bach1", "Ifi204", "Oasl2", "Ifi203", "Parp14", "Gbp3",
    "Slfn2", "Ptprc", "Cyrib", "Tgfbi", "Gpr141", "Ifit2", "Cd300lf", "Gbp7", "Rnf213", "Fcgr1",
    "Fgl2", "Atp6v0c", "Ifi205", "Gbp4", "Cebpb", "Lilrb4b", "Tgm2", "Pkm", "Gbp2", "Ifi211",
    "Lilrb4a", "Cxcl10", "Ctss", "Ly6a"
  ),
  CrossPresentation_DC_Strict = c(
    "Akap13", "Arl4c", "Fam118a", "Fos", "Fosb", "H2-Q7", "H2-Q6", "H2-Q4", "H2-Q10", "H2-Q2",
    "H2-K1", "H2-T23", "H2-D1", "H2-Aa", "H2-Ab1", "Icam1", "Irf1", "Irf7", "Isg20", "Jun",
    "Junb", "Klf6", "Ldlrad4", "Lilra6", "Pira12", "Mx1", "Pde4b", "Plac8", "Plek", "Ranbp2",
    "Sell", "Snhg5", "Stat1", "Tcl1", "Tspyl2", "Ucp2", "Wars1"
  ),
  CrossPresentation_Macrophage_Strict = c(
    "Apoc1", "Apoe", "B2m", "C1qa", "C1qb", "C1qc", "C3", "Ccl3", "Ccl4", "Cxcl10", "Cxcl9",
    "Fcgr4", "Fn1", "H2-Q7", "H2-Q6", "H2-Q4", "H2-Q10", "H2-Q2", "H2-K1", "H2-T23", "H2-D1",
    "H2-Aa", "H2-Eb1", "Hspa8", "Ier2", "Ifi27", "Ifitm1", "Ifitm3", "Lgals3bp", "Ly6e", "Ncf1",
    "Psap", "Psme2", "Serping1", "Stat1", "Vamp5", "Wars1"
  )
)

human_signatures <- list(
  MouseToHuman_ISG47 = c(
    "ARID5B", "ATP11B", "ATP6V0C", "BACH1", "BATF2", "CD300LF", "CEBPB", "CREB5",
    "CTSS", "CXCL10", "CYRIB", "ENTPD1", "FCGR1A", "FCGR3B", "FCGR3A", "FGL2",
    "GBP2", "GNA13", "GPR141", "PYHIN1", "MNDA", "IFI27L2", "IFIT2", "IRGM",
    "IL10RA", "IRAK2", "IRF7", "LILRB4", "OAS3", "PARP14", "PDE7B", "PIK3AP1",
    "PNP", "PTPRC", "RNF213", "RTP4", "SAMD9L", "SAMHD1", "SDCBP", "SLAMF8",
    "SLFN13", "STAT2", "TGFBI", "TGM2", "TNFAIP2", "XAF1", "ZBP1"
  )
)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1") && !file.exists(mouse_rds) && !file.exists(human_rds)) {
  message("SMOKE_TEST=1: no mouse or human Seurat object found; skipping UCell confirmation.")
  quit(save = "no", status = 0)
}

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: UCell version ", as.character(packageVersion("UCell")))
  message("SMOKE_TEST=1: mouse signatures defined: ", length(mouse_signatures))
  message("SMOKE_TEST=1: human signatures defined: ", length(human_signatures))
  message("SMOKE_TEST=1: skipping full UCell scoring.")
  quit(save = "no", status = 0)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

layer_matrix <- function(object, assay = "RNA") {
  DefaultAssay(object) <- assay
  layer <- if (ncol(GetAssayData(object, assay = assay, slot = "counts")) > 0) "counts" else "data"
  list(
    matrix = GetAssayData(object, assay = assay, slot = layer),
    layer = layer
  )
}

gene_presence <- function(signatures, features, species) {
  bind_rows(lapply(names(signatures), function(signature) {
    genes <- signatures[[signature]]
    data.frame(
      species = species,
      signature = signature,
      gene = genes,
      status = ifelse(genes %in% features, "Present", "Missing"),
      stringsAsFactors = FALSE
    )
  }))
}

score_object <- function(object, signatures, species, outdir) {
  if (!assay_name %in% Assays(object)) stop("Missing assay in ", species, " object: ", assay_name)

  expr <- layer_matrix(object, assay = assay_name)
  present_signatures <- lapply(signatures, intersect, rownames(expr$matrix))
  present_signatures <- present_signatures[lengths(present_signatures) > 0]
  if (length(present_signatures) == 0) stop("No signature genes are present in ", species, " object.")

  scores <- as.data.frame(UCell::ScoreSignatures_UCell(
    matrix = expr$matrix,
    features = present_signatures,
    name = "_UCell",
    ncores = as.integer(Sys.getenv("UCELL_NCORES", unset = "1"))
  ))
  scores$cell <- rownames(scores)

  metadata <- object@meta.data %>%
    mutate(cell = rownames(object@meta.data))

  cell_scores <- metadata %>%
    left_join(scores, by = "cell")

  score_cols <- setdiff(colnames(scores), "cell")
  summary_vars <- intersect(c(cluster_col, sample_col, condition_col, "seurat_clusters", "patient", "response", "timepoint"), colnames(cell_scores))

  cell_summary <- cell_scores %>%
    pivot_longer(all_of(score_cols), names_to = "signature", values_to = "ucell_score") %>%
    group_by(across(all_of(c("signature", summary_vars)))) %>%
    summarise(
      n_cells = n(),
      mean_score = mean(ucell_score, na.rm = TRUE),
      median_score = median(ucell_score, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(cell_scores, file.path(outdir, paste0(species, "_ucell_cell_scores.csv")), row.names = FALSE)
  write.csv(cell_summary, file.path(outdir, paste0(species, "_ucell_score_summary.csv")), row.names = FALSE)

  if (cluster_col %in% colnames(cell_scores) && condition_col %in% colnames(cell_scores)) {
    plot_df <- cell_scores %>%
      pivot_longer(all_of(score_cols), names_to = "signature", values_to = "ucell_score") %>%
      mutate(cluster = as.character(.data[[cluster_col]]))

    p <- ggplot(plot_df, aes(x = cluster, y = ucell_score, fill = .data[[condition_col]])) +
      geom_boxplot(outlier.shape = NA, width = 0.65) +
      facet_wrap(~signature, scales = "free_y") +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(x = "Cluster", y = "UCell score", fill = "Condition")

    ggsave(file.path(outdir, paste0(species, "_ucell_scores_by_cluster_condition.pdf")), p, width = 11, height = 7)
  }

  list(scores = cell_scores, summary = cell_summary, layer = expr$layer)
}

tables <- list()

if (file.exists(mouse_rds)) {
  mouse_obj <- readRDS(mouse_rds)
  tables$mouse_gene_presence <- gene_presence(mouse_signatures, rownames(mouse_obj[[assay_name]]), "mouse")
  mouse_results <- score_object(mouse_obj, mouse_signatures, "mouse", output_dir)
  tables$mouse_summary <- mouse_results$summary
}

if (file.exists(human_rds)) {
  human_obj <- readRDS(human_rds)
  tables$human_gene_presence <- gene_presence(human_signatures, rownames(human_obj[[assay_name]]), "human")
  human_results <- score_object(human_obj, human_signatures, "human", output_dir)
  tables$human_summary <- human_results$summary
}

if (length(tables) == 0) stop("No input Seurat objects found.")

write.xlsx(tables, file.path(output_dir, "ucell_signature_confirmation.xlsx"), overwrite = TRUE)
message("Done. UCell confirmation outputs written to: ", output_dir)
