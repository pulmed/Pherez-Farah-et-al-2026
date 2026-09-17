#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/cross_species/translate_mouse_signatures_to_human.R
# Original file: translate_genes.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Translate mouse signature genes to human orthologs for cross-species interpretation.
# Inputs: Final mouse signature gene vectors defined in this script.
# Outputs: Mouse-to-human ortholog CSV.
# Assay/layer input: Not applicable.
# Dependencies: dplyr, tidyr, tibble, orthogene.
# Environment: Main analysis environment (conda + renv), plus optional orthogene/gprofiler support.
# Notes:
# - This helper uses final repository signature definitions where updated versions exist.
# - The orthogene/gprofiler conversion may require an internet connection.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

source("scripts/utils/seurat_io.R")

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- make_output_dir(file.path(output_root, "cross_species", "ortholog_translation"))
output_csv <- Sys.getenv("ORTHOLOG_OUTPUT_CSV", unset = file.path(outdir, "mouse_to_human_orthologues.csv"))
ortholog_method <- Sys.getenv("ORTHOLOG_METHOD", unset = "gprofiler")
non121_strategy <- Sys.getenv("ORTHOLOG_NON121_STRATEGY", unset = "keep_both_species")

# ------------------------------------------------------------------------------
# FINAL/REPOSITORY SIGNATURE DEFINITIONS
# ------------------------------------------------------------------------------
isg_genes <- c(
  "Irgm2", "Fcgr3", "Il10ra", "Pik3ap1", "Oas3", "Ifi47", "Entpd1", "Tnfaip2", "Rtp4", "Gna13",
  "Irgm1", "Fcgr4", "Stat2", "Sdcbp", "Gm4951", "Ifi213", "Samd9l", "Igtp", "Ifi207", "Xaf1",
  "Socs3", "Mndal", "Pde7b", "Irf7", "Gm12185", "Slamf8", "Atp11b", "Samhd1", "Fbxl5", "Tgtp2",
  "Slfn1", "Cd86", "Ifi27l2a", "Gm20663", "Zbp1", "Slfn8", "Creb5", "Irak2", "Ifi206", "Pnp",
  "Iigp1", "Ly6i", "Batf2", "Arid5b", "Bach1", "Ifi204", "Oasl2", "Ifi203", "Parp14", "Gbp3",
  "Slfn2", "Ptprc", "Cyrib", "Tgfbi", "Gpr141", "Ifit2", "Cd300lf", "Gbp7", "Rnf213", "Fcgr1",
  "Fgl2", "Atp6v0c", "Ifi205", "Gbp4", "Cebpb", "Lilrb4b", "Tgm2", "Pkm", "Gbp2", "Ifi211",
  "Lilrb4a", "Cxcl10", "Ctss", "Ly6a"
)

m1_genes <- c("Nos2", "Il1a", "Il1b", "Tnf", "Il6", "Cd80", "Cd86", "Cxcl9", "Cxcl10", "Stat1")
m2_genes <- c("Arg1", "Mrc1", "Retnla", "Chil3", "Il10", "Tgfb1", "Cd163", "Ccl17", "Ccl22", "Stat6")

cd301b_dcs <- c(
  "Pxdc1", "Ass1", "Vdr", "Pmvk", "Arhgef10", "Cish", "Upp1", "H2-Oa", "Lsr",
  "Mgl2", "Ankrd37", "H2-DMb2", "Slc4a11", "Il1r2", "Rhoq", "Klrb1b", "Fgl2",
  "Klrk1", "Afdn", "Ifitm1", "Dapk1", "Klrd1", "Adgrg5", "Itgax", "Spint1"
)

lam_genes <- c("Trem2", "Lipa", "Lpl", "Ctsb", "Ctsl", "Fabp4", "Fabp5", "Lgals1", "Lgals3", "Cd9", "Cd36")

cross_presentation_genes <- c(
  "H2-K1", "H2-D1", "B2m",
  "Tap1", "Tap2", "Tapbp",
  "Calr", "Canx", "Pdia3", "Erap1",
  "Psmb8", "Psmb9", "Psmb10", "Psme1", "Psme2",
  "Sec61a1", "Sec61b", "Sec61g",
  "Cybb", "Rab27a", "Rac2", "Vamp8", "Stx4",
  "Wdfy4", "Sec22b", "Lnpep"
)

signatures <- list(
  ISG = isg_genes,
  M1_final = m1_genes,
  M2_final = m2_genes,
  CD301b_DCs = cd301b_dcs,
  LAM = lam_genes,
  CrossPresentation_final = cross_presentation_genes
)

input_table <- enframe(signatures, name = "signature", value = "mouse_gene") %>%
  unnest(mouse_gene) %>%
  distinct(signature, mouse_gene)

if (identical(Sys.getenv("SMOKE_TEST", unset = "0"), "1")) {
  message("SMOKE_TEST=1: signatures: ", length(signatures))
  message("SMOKE_TEST=1: unique mouse genes to translate: ", length(unique(input_table$mouse_gene)))
  quit(save = "no", status = 0)
}

if (!requireNamespace("orthogene", quietly = TRUE)) {
  stop(
    "Package 'orthogene' is required for ortholog translation but is not installed. ",
    "Install it in the analysis environment before running this optional helper.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# ORTHOLOG CONVERSION
# ------------------------------------------------------------------------------
all_mouse_genes <- unique(input_table$mouse_gene)

map_raw <- orthogene::convert_orthologs(
  gene_df = data.frame(mouse_gene = all_mouse_genes),
  gene_input = "mouse_gene",
  gene_output = "columns",
  input_species = "mouse",
  output_species = "human",
  method = ortholog_method,
  non121_strategy = non121_strategy
)

required_columns <- c("input_gene", "ortholog_gene")
if (!all(required_columns %in% colnames(map_raw))) {
  stop(
    "Unexpected orthogene output columns: ",
    paste(colnames(map_raw), collapse = ", "),
    call. = FALSE
  )
}

map_table <- map_raw %>%
  select(mouse_gene = input_gene, human_gene = ortholog_gene) %>%
  distinct()

final_table <- input_table %>%
  left_join(map_table, by = "mouse_gene") %>%
  mutate(
    status = if_else(is.na(human_gene), "not_converted", "orthogene"),
    method = ortholog_method,
    non121_strategy = non121_strategy
  ) %>%
  arrange(signature, mouse_gene, human_gene)

write.csv(final_table, output_csv, row.names = FALSE)

summary_table <- final_table %>%
  group_by(signature) %>%
  summarise(
    n_mouse_genes = n_distinct(mouse_gene),
    n_mouse_genes_with_human_ortholog = n_distinct(mouse_gene[!is.na(human_gene)]),
    n_human_ortholog_rows = sum(!is.na(human_gene)),
    .groups = "drop"
  )

write.csv(summary_table, file.path(outdir, "mouse_to_human_orthologue_summary.csv"), row.names = FALSE)

message("Done. Ortholog table written to: ", output_csv)
