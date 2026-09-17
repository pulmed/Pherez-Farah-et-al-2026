#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/confirm_signature_trends_sample_level.R
# Original file: 20260330 SAMPLE LEVEL CONFIRMATION OF SIGNATURE TRENDS.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Confirm selected signature trends at sample level using heatmaps and ridgeplots.
# Inputs: Canonical analysis-ready Seurat object with RNA data, final_clusters, and hash.ID metadata.
# Outputs: Sample-level signature heatmaps, ridgeplots, summary CSVs, and text summaries.
# Assay/layer input: RNA data.
# Dependencies: Seurat, dplyr, tidyr, ggplot2, ggridges, pheatmap.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - The script contains repeated signature-specific blocks preserved from the original analysis.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "MREG"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("13","1","7","8","9")

# ------------------------------------------------------------------------------
# FULL MREG SIGNATURE
# ------------------------------------------------------------------------------
sig_genes <- c(
  "Slco5a1","Ly96","Ogfrl1","Kansl3","Tmem131","Tbc1d8","Map4k4","Nabp1","Stat4","Stat1","Cflar","Icos","Nrp2","Ino80d",
  "Gm20342","Adam23","Idh1","Mreg","Dnajb2","C130026I21Rik","A530032D15Rik","Sp140","Sh3bp4","Agap1","Gpc1","Tnfrsf11a",
  "Phlpp1","Nuak2","Gm38399","Csrp1","Phlda3","Lad1","Cacna1s","Kif21b","Rgs1","Cep350","Rabgap1l","Rcsd1","Pogk","Uap1",
  "Slamf1","Pyhin1","Mfsd7b","G0s2","Dclre1c","Frmd4a","Pfkfb3","Gm10851","Il2ra","Il15ra","Il1rn","Traf2","Tsc1","Ass1",
  "St6galnac6","Stxbp1","Traf1","Ggta1","Arl5a","Cytip","Gm13546","Ly75","Tank","Gca","Cers6","Nostrin","Itga4","Dnajc10",
  "Ube2l6","Arhgap1","Phf21a","Traf6","D330050G23Rik","Spred1","Bahd1","Pla2g4f","Vps39","B2m","Sema6d","Hdc","Tmem127",
  "Stard7","Prnp","Rassf2","Gpcpd1","Jag1","Flrt3","Dstn","Ralgapa2","Srxn1","Procr","Src","Stk4","Slpi","Cd40","Zmynd8",
  "Rab22a","Ppdpf","Helz2","Arfrp1","Kdm6a","Chst7","Enox2","Renbp","Emd","Gdi1","Apool","Rab9","Gm16685","Fabp5","Fabp4",
  "Car13","Car2","Gyg","Gnb4","Slc7a11","Rap2b","E130311K13Rik","Slc33a1","1110032F04Rik","Etfdh","D930015E06Rik","Cd1d1",
  "Etv3","Pmvk","Slc27a3","Ints3","Npr1","S100a8","S100a9","Golph3l","Ankrd35","Pias3","Nudt17","Polr3c","Rnf115","Pde4dip",
  "Igsf3","Mab21l3","Rhoc","Gnai3","Stxbp3","Vcam1","Plppr4","Fnbp1l","Synpo2","Dkk2","Gbp5","Gbp2","Pnrc1","3110043O21Rik",
  "Glipr2","Nans","Galnt12","Nr4a3","Invs","Akap2","Col27a1","Cdkn2a","Cdkn2b","Nsun4","St3gal3","9530034E10Rik","Kdm4a"
)

# ------------------------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
stopifnot(sample_col %in% colnames(obj[[]]))

obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

# ------------------------------------------------------------------------------
# FORCE NORMALIZATION
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

present_genes <- intersect(sig_genes, rownames(expr))
missing_genes <- setdiff(sig_genes, present_genes)

write.csv(
  data.frame(gene = present_genes),
  file.path(outdir, "MREG_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing_genes),
  file.path(outdir, "MREG_missing_genes.csv"),
  row.names = FALSE
)

if (length(present_genes) == 0) {
  stop("No MREG genes found after normalization.")
}

# ------------------------------------------------------------------------------
# FILTER CELLS
# ------------------------------------------------------------------------------
md <- obj[[]]

cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    md[[cluster_col]] %in% clusters_use
]

obj <- subset(obj, cells = cells_keep)
md <- obj[[]]
expr <- expr[present_genes, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
scale_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  z
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

mean_by_group <- function(mat, groups, ordered_levels) {
  out <- lapply(ordered_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(mat)
  colnames(out) <- ordered_levels
  out
}

# ------------------------------------------------------------------------------
# SAMPLE LOOP
# ------------------------------------------------------------------------------
samples_use <- sort(unique(md[[sample_col]]))

for (s in samples_use) {

  message("Processing sample: ", s)

  cells_s <- rownames(md)[md[[sample_col]] == s]
  if (length(cells_s) == 0) next

  md_s <- md[cells_s, , drop = FALSE]
  clusters_present <- clusters_use[clusters_use %in% unique(md_s[[cluster_col]])]

  if (length(clusters_present) == 0) next

  expr_s <- expr[, cells_s, drop = FALSE]

  avg_expr <- mean_by_group(
    mat = expr_s,
    groups = md_s[[cluster_col]],
    ordered_levels = clusters_present
  )

  zmat <- scale_rows(avg_expr)
  sample_tag <- safe_name(s)

  write.csv(
    as.data.frame(avg_expr),
    file.path(outdir, paste0(sample_tag, "_avg_expression.csv"))
  )

  write.csv(
    as.data.frame(zmat),
    file.path(outdir, paste0(sample_tag, "_row_zscore.csv"))
  )

  png(
    file.path(outdir, paste0(sample_tag, "_MREG_heatmap.png")),
    width = 2400,
    height = max(2200, 10 * nrow(zmat)),
    res = 300
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 4,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("MREG | ", s)
  )

  dev.off()

  pdf(
    file.path(outdir, paste0(sample_tag, "_MREG_heatmap.pdf")),
    width = 8,
    height = max(10, nrow(zmat) * 0.08)
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 4,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("MREG | ", s)
  )

  dev.off()
}

message("DONE: all sample-level MREG heatmaps saved in:")
message(outdir)

# confirmation of ridge

# RIDGEPLOT: MREG signature across clusters, SAMPLE as replication unit
# Same concept as original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/MREG_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "MREG_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("13","1","7","8","9")

sig_genes <- c(
  "Slco5a1","Ly96","Ogfrl1","Kansl3","Tmem131","Tbc1d8","Map4k4","Nabp1","Stat4","Stat1","Cflar","Icos","Nrp2","Ino80d",
  "Gm20342","Adam23","Idh1","Mreg","Dnajb2","C130026I21Rik","A530032D15Rik","Sp140","Sh3bp4","Agap1","Gpc1","Tnfrsf11a",
  "Phlpp1","Nuak2","Gm38399","Csrp1","Phlda3","Lad1","Cacna1s","Kif21b","Rgs1","Cep350","Rabgap1l","Rcsd1","Pogk","Uap1",
  "Slamf1","Pyhin1","Mfsd7b","G0s2","Dclre1c","Frmd4a","Pfkfb3","Gm10851","Il2ra","Il15ra","Il1rn","Traf2","Tsc1","Ass1",
  "St6galnac6","Stxbp1","Traf1","Ggta1","Arl5a","Cytip","Gm13546","Ly75","Tank","Gca","Cers6","Nostrin","Itga4","Dnajc10",
  "Ube2l6","Arhgap1","Phf21a","Traf6","D330050G23Rik","Spred1","Bahd1","Pla2g4f","Vps39","B2m","Sema6d","Hdc","Tmem127",
  "Stard7","Prnp","Rassf2","Gpcpd1","Jag1","Flrt3","Dstn","Ralgapa2","Srxn1","Procr","Src","Stk4","Slpi","Cd40","Zmynd8",
  "Rab22a","Ppdpf","Helz2","Arfrp1","Kdm6a","Chst7","Enox2","Renbp","Emd","Gdi1","Apool","Rab9","Gm16685","Fabp5","Fabp4",
  "Car13","Car2","Gyg","Gnb4","Slc7a11","Rap2b","E130311K13Rik","Slc33a1","1110032F04Rik","Etfdh","D930015E06Rik","Cd1d1",
  "Etv3","Pmvk","Slc27a3","Ints3","Npr1","S100a8","S100a9","Golph3l","Ankrd35","Pias3","Nudt17","Polr3c","Rnf115","Pde4dip",
  "Igsf3","Mab21l3","Rhoc","Gnai3","Stxbp3","Vcam1","Plppr4","Fnbp1l","Synpo2","Dkk2","Gbp5","Gbp2","Pnrc1","3110043O21Rik",
  "Glipr2","Nans","Galnt12","Nr4a3","Invs","Akap2","Col27a1","Cdkn2a","Cdkn2b","Nsun4","St3gal3","9530034E10Rik","Kdm4a",
  "Inpp5b","Mtf1","Gnl2","Ago3","Ago1","Rnf19b","Bsdc1","Marcksl1","Epb41","Zdhhc18","Extl1","Clic4","Sh2d5","Iffo2",
  "Plekhm2","Vps13d","Tnfrsf1b","Kif1b","Spsb1","Tnfrsf9","Mmp23","Fam132a","Tnfrsf4","Isg15","Krit1","Sri","Nub1","Tmem128",
  "Sepsecs","Pi4k2b","Arap2","Tbc1d1","Gm42726","Nipal1","Rufy3","Cxcl1","Anxa3","Bmp2k","Tmem150c","Klhl8","Nudt9","Gbp8",
  "Gbp9","Gbp4","Pxmp2","Adrbk2","Iscu","Ssh1","Oasl1","Gm10399","Vsig10","Tesc","Hvcn1","Tmem120b","Sbno1","Asl","Rabgef1",
  "Gatsl2","Ncf1","Clip2","Wbscr27","Rasa4","Adap1","Fbxl18","Fscn1","Cyth3","Smurf1","Lnx2","N4bp2l1","Tmem168",
  "B630005N14Rik","Tes","Ahcyl2","Strip2","Tmem140","Zfp467","Herc6","Vopp1","Ndnf","Il12rb2","Cd8a","M1ap","Dok1","Loxl3",
  "Htra2","Mthfd2","Spr","Asprv1","Mxd1","Anxa4","Gfpt1","Slc6a6","Ppp4r2","Wnk1","Mical3","Usp18","Apobec1","Eno2","Tapbpl",
  "Clec2i","Clec2d","Bcl2l14","Crebl2","Hebp1","Aebp2","Etnk1","Caprin2","Fam60a","3010003L21Rik","Amn1","Zfp524","Ceacam15",
  "Zfp296","Relb","Clptm1","Ceacam1","Sertad1","Plekhg2","Nfkbib","Kcnk6","2200002D01Rik","Spint2","Alkbh6","Sdhaf1","Hspb6",
  "Nkg7","Spib","Atf5","Il4i1","Ftl1","Rasip1","Ndnl2","Lrrk1","RP23-49I18.5","Arpin","Mex3b","Zfand6","Prkrir","Stard10",
  "Il18bp","Pgap2","Swap70","Sbf2","Rras2","Gga2","Il21r","Mvp","Cd2bp2","Prr14","Fam53b","Adam8","Ptdss2","Stx11","Adgrg6",
  "Ncoa7","Marcks","Hsf2","Sgpl1","Gm5424","Dnajc12","Nrbf2","Ccdc6","Adora2a","Gucd1","Gstt1","Icosl","Cstb","Cnn2","Efna2",
  "Mob3a","Izumo4","Gadd45b","Syn3","Dram1","Chpt1","Cdk17","Tmcc3","Plxnc1","Socs2","Btg1","Tmtc2","Osbpl8","Tbc1d15","Rab21",
  "Tmem19","Thap2","Zfc3h1","Tmbim4","Rassf3","Ddit3","Ikzf4","Suox","Cd63","Rab20","Atp11a","Lamp1","Polb","Nrg1","Ubxn8",
  "Gsr","Mfhas1","Casp3","Nr2f6","Fam129c","Fam32a","Large","Il15","Hook2","9330175E14Rik","Nlrc5","Cpne2","Ccl22","Cx3cl1",
  "Gins3","Ndrg4","Nae1","Rrad","Nfat5","Hp","Mon1b","Cmc2","Atmin","Crispld2","Fam89a","Itgb1","Arf4","Arhgap22","Ccser2",
  "Ghitm","Ddhd1","Ktn1","Ndrg2","Arhgef40","Psme2","Kpna3","Clu","Epsti1","Sugt1","Tbc1d4","Tmtc4","Tmem123","Birc2","Birc3",
  "Fut4","Icam1","Zfp809","Vwa5a","Gramd1b","Ddx6","Phldb1","Scn2b","Fxyd2","Dscaml1","Sik3","Zc3h12c","Tspan3","Peak1","Mpi",
  "Sema7a","Stoml1","Dennd4a","Rab8b","Rps27l","Lactb","Tpm1","Rora","Aqp9","Aldh1a2","Gclc","Hmgn3","Mthfsl","Bcl2a1d","Bcl2a1a",
  "Bcl2a1b","Rasa2","Nmnat3","Pcbp4","Uba7","Bcl2a1c","Wdr48","9530059O14Rik","Tmem158","Gatsl3","Myo1g","Ramp3","Aftph","Peli1",
  "Pex13","Pus10","Rel","Papolg","Ccdc88a","Cpeb4","Il12b","Trim7","Psme2b","Sqstm1","Ltc4s","Tcf7","Irf1","Gm12216","Pdlim4",
  "Specc1","Pik3r5","Dnah2","Cxcl16","Zmynd15","Tm4sf5","Eno3","Nup88","Txndc17","Ccl5","Mtmr4","Bzrap1","Mmd","Gngt2","Nfe2l1",
  "Mllt6","Arl5c","Ccr7","Smarce1","Stat5a","Stat3","Adam11","Map3k14","Ace","Nup85","Gga3","Mif4gd","Mrpl38","Cyth1","Rptor",
  "Slc25a10","Net1","Klf6","Idi1","Gtpbp4","Gpr137b-ps","Gpr137b","Sox4","Cdkal1","Serpinb1a","Serpinb6b","Serpinb9","Serpinb6a",
  "Tubb2b","Psmg4","Slc22a23","Hivep1","Cd83","Mylip","Atxn1","Kdm1b","Auh","Nfil3","Cdc14b","Gm10116","Serinc5","Pik3r1","Pde4d",
  "Map3k1","Gzma","Fam49a","Nampt","Gdap10","Atxn7l1","Scin","Pnpla8","Fam177a","Sav1","Daam1","Plek2","Zfp36l1","2310015A10Rik",
  "Ttll5","Ift43","Gtf2a1","Nrde2","Calm1","Slc25a29","Traf3","Exoc3l4","Inf2","Gpr132","Itgb8","Capsl","Il7r","Basp1","Ankrd33b",
  "Laptm4b","Rnf19a","Ndrg1","Arc","Ppp1r16a","Apol7c","Apol10b","C1qtnf6","Triobp","Micall1","Josd1","Gtpbp1","Apobec3","Zc3h7b",
  "Slc38a2","Adcy6","Cacnb3","Rhebl1","Tuba1a","Prph","Lima1","Dip2b","Slc4a8","Grasp","Rogdi","Socs1","Snn","Dnm1l","B3gnt5",
  "Yeats2","Ccdc50","Pcyt1a","Rubcn","Sec22a","Hspbap1","Pla1a","Cd80","Poglut1","Tmem39a","Arhgap31","Cd200","Retnlg","Ift57","Cblb",
  "Filip1l","Samsn1","Nrip1","Mir155hg","Mx1","Rnaset2a","Mmp25","Rpl3l","Tmem8","D17Wsu92e","H2-K1","H2-Eb2","H2-Q4","H2-Q6","H2-Q7",
  "H2-M2","Tmem63b","Mrpl14","Gtpbp2","Foxp4","Rftn1","Stap2","Cd70","Trip10","Wash1","Arhgap28","Cdc42ep3","Fbxo11","Foxn2","Crem","Cdh2",
  "Gypc","Alpk2","Chka","Cdc42ep2","Ehd1","Atl3","Cpsf7","AW112010","Ms4a7","Nmrk1","Jak2","Insl6","Plgrkt","Cd274","Pdcd1lg2","Ric1","Papss2",
  "Fas","Pcgf5","Fgfbp3","March5","Avpi1","Got1","Nfkb2","Gsto1","4833407H14Rik","Dusp5","Vti1a","Tcf7l2","Dclre1a","Atrnl1","PISD","Umad1"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "MREG_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(plot_df, file.path(outdir, "MREG_sample_level_scores_by_cluster.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(plot_df, aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "MREG signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(file.path(outdir, "MREG_sample_replicate_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "MREG_sample_replicate_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "MREG_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)

# SPRANGER

# FINAL CLEAN SAMPLE-LEVEL HEATMAPS: SPRANGER SIGNATURE
# Robust version
# - forces RNA log-normalization
# - one heatmap per sample (hash.ID)
# - averages normalized expression by cluster inside each sample
# - row z-score only for visualization
# - Seurat v5 safe
# ------------------------------------------------------------------------------
# INPUT:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# OUTPUT:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/SPRANGER

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "SPRANGER"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("13","1","7","8","9")

# ------------------------------------------------------------------------------
# SPRANGER SIGNATURE
# ------------------------------------------------------------------------------
sig_genes <- c(
  "Ifit1bl1","Ifit3b","Ifit3","Ifit1","Cxcl10","Ifit2","Rsad2","Phf11d","Isg20","Gbp7","Slfn1","Ifi213","Fam26f","Ifi47",
  "Ifi205","Nt5c3","Gbp2","Ifih1","Ifi206","Oasl1","Ifi211","Fgl2","Ddx58","1600014C10Rik","Slfn9","Igtp","Cmpk2","Lipg",
  "Herc6","Ms4a6d","Ms4a4b","Endod1","Slfn8","Ifi214","Ddx60","Dhx58","Tor3a","Iigp1","A530064D06Rik","Znfx1","Axl","Ifi44",
  "Mx1","Oas3","Gbp3","Carhsp1","Fcgr4","Gbp4","Slfn4","Treml2","Mthfr","Lpxn","Adap2","Phf11c","Svbp","AA467197","Trim30b",
  "B430306N03Rik","Cd69","Rnpep","Ifi204","Ms4a4c","Isg15","Fcgr1","Irf7","Phf11b","Rtp4","Usp18","Ly6a","Ms4a6b","Ifi209",
  "Mndal","Xaf1","Ccnd2","Oasl2","Zbp1","Slfn5","Ifi203","Ly6e","Samd9l","Phf11a","Ms4a6c","Slfn2","Sp100","H2-T22","Samhd1",
  "Ctss","Ifitm3","Sat1","Ube2l6","Lgals3bp","Parp14","Pnp","Sp110","Selenow","Sdcbp","Snx2","Irgm1","Tspo","Trafd1","Stat2",
  "Daxx","Ifi207","Fcer1g","Hck","Dck","Pttg1","Chmp4b","Rnf213","Nmi","Klrk1","Hmox2","Aftph","Lgals9","Ifi27l2a","Rnf34",
  "Usp25","Tmem219","Wdfy1","2810474O19Rik","Ccr5","Scimp","Pstpip1","Dtx3l","Eif2ak2","Lilr4b","Ctsc","Cd86"
)

# ------------------------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
stopifnot(sample_col %in% colnames(obj[[]]))

obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

# ------------------------------------------------------------------------------
# FORCE NORMALIZATION
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

present_genes <- intersect(sig_genes, rownames(expr))
missing_genes <- setdiff(sig_genes, present_genes)

write.csv(
  data.frame(gene = present_genes),
  file.path(outdir, "SPRANGER_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing_genes),
  file.path(outdir, "SPRANGER_missing_genes.csv"),
  row.names = FALSE
)

if (length(present_genes) == 0) {
  stop("No SPRANGER genes found after normalization.")
}

# ------------------------------------------------------------------------------
# FILTER CELLS
# ------------------------------------------------------------------------------
md <- obj[[]]

cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    md[[cluster_col]] %in% clusters_use
]

obj <- subset(obj, cells = cells_keep)
md <- obj[[]]
expr <- expr[present_genes, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
scale_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  z
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

mean_by_group <- function(mat, groups, ordered_levels) {
  out <- lapply(ordered_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(mat)
  colnames(out) <- ordered_levels
  out
}

# ------------------------------------------------------------------------------
# SAMPLE LOOP
# ------------------------------------------------------------------------------
samples_use <- sort(unique(md[[sample_col]]))

for (s in samples_use) {

  message("Processing sample: ", s)

  cells_s <- rownames(md)[md[[sample_col]] == s]
  if (length(cells_s) == 0) next

  md_s <- md[cells_s, , drop = FALSE]
  clusters_present <- clusters_use[clusters_use %in% unique(md_s[[cluster_col]])]

  if (length(clusters_present) == 0) next

  expr_s <- expr[, cells_s, drop = FALSE]

  avg_expr <- mean_by_group(
    mat = expr_s,
    groups = md_s[[cluster_col]],
    ordered_levels = clusters_present
  )

  zmat <- scale_rows(avg_expr)
  sample_tag <- safe_name(s)

  write.csv(
    as.data.frame(avg_expr),
    file.path(outdir, paste0(sample_tag, "_avg_expression.csv"))
  )

  write.csv(
    as.data.frame(zmat),
    file.path(outdir, paste0(sample_tag, "_row_zscore.csv"))
  )

  png(
    file.path(outdir, paste0(sample_tag, "_SPRANGER_heatmap.png")),
    width = 2400,
    height = max(2200, 10 * nrow(zmat)),
    res = 300
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("SPRANGER | ", s)
  )

  dev.off()

  pdf(
    file.path(outdir, paste0(sample_tag, "_SPRANGER_heatmap.pdf")),
    width = 8,
    height = max(10, nrow(zmat) * 0.14)
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("SPRANGER | ", s)
  )

  dev.off()
}

message("DONE: all sample-level SPRANGER heatmaps saved in:")
message(outdir)

# Confirmation of rige, spranger

# RIDGEPLOT: SPRANGER signature across clusters, SAMPLE as replication unit
# Same concept as original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/SPRANGER_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "SPRANGER_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("13","1","7","8","9")

sig_genes <- c(
  "Ifit1bl1","Ifit3b","Ifit3","Ifit1","Cxcl10","Ifit2","Rsad2","Phf11d","Isg20","Gbp7","Slfn1","Ifi213","Fam26f","Ifi47","Ifi205","Nt5c3",
  "Gbp2","Ifih1","Ifi206","Oasl1","Ifi211","Fgl2","Ddx58","1600014C10Rik","Slfn9","Igtp","Cmpk2","Lipg","Herc6","Ms4a6d","Ms4a4b",
  "Endod1","Slfn8","Ifi214","Ddx60","Dhx58","Tor3a","Iigp1","A530064D06Rik","Znfx1","Axl","Ifi44","Mx1","Oas3","Gbp3","Carhsp1",
  "Fcgr4","Gbp4","Slfn4","Treml2","Mthfr","Lpxn","Adap2","Phf11c","Svbp","AA467197","Trim30b","B430306N03Rik","Cd69","Rnpep","Ifi204",
  "Ms4a4c","Isg15","Fcgr1","Irf7","Phf11b","Rtp4","Usp18","Ly6a","Ms4a6b","Ifi209","Mndal","Xaf1","Ccnd2","Oasl2","Zbp1","Slfn5",
  "Ifi203","Ly6e","Samd9l","Phf11a","Ms4a6c","Slfn2","Sp100","H2-T22","Samhd1","Ctss","Ifitm3","Sat1","Ube2l6","Lgals3bp","Parp14",
  "Pnp","Sp110","Selenow","Sdcbp","Snx2","Irgm1","Tspo","Trafd1","Stat2","Daxx","Ifi207","Fcer1g","Hck","Dck","Pttg1","Chmp4b","Rnf213",
  "Nmi","Klrk1","Hmox2","Aftph","Lgals9","Ifi27l2a","Rnf34","Usp25","Tmem219","Wdfy1","2810474O19Rik","Ccr5","Scimp","Pstpip1",
  "Dtx3l","Eif2ak2","Lilr4b","Ctsc","Cd86"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "SPRANGER_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(plot_df, file.path(outdir, "SPRANGER_sample_level_scores_by_cluster.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(plot_df, aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "SPRANGER signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(file.path(outdir, "SPRANGER_sample_replicate_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "SPRANGER_sample_replicate_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "SPRANGER_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)

# LAM heatmaps per sample

# FINAL CLEAN SAMPLE-LEVEL HEATMAPS: LAM SIGNATURE
# Robust version
# - forces RNA log-normalization
# - one heatmap per sample (hash.ID)
# - averages normalized expression by cluster inside each sample
# - row z-score only for visualization
# - Seurat v5 safe
# ------------------------------------------------------------------------------
# INPUT:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# OUTPUT:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/LAM

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "LAM"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14")

# ------------------------------------------------------------------------------
# LAM SIGNATURE
# ------------------------------------------------------------------------------
sig_genes <- c(
  "Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5","Lgals1","Lgals3","Cd9","Cd36"
)

# ------------------------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
stopifnot(sample_col %in% colnames(obj[[]]))

obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

# ------------------------------------------------------------------------------
# FORCE NORMALIZATION
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

present_genes <- intersect(sig_genes, rownames(expr))
missing_genes <- setdiff(sig_genes, present_genes)

write.csv(
  data.frame(gene = present_genes),
  file.path(outdir, "LAM_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing_genes),
  file.path(outdir, "LAM_missing_genes.csv"),
  row.names = FALSE
)

if (length(present_genes) == 0) {
  stop("No LAM genes found after normalization.")
}

# ------------------------------------------------------------------------------
# FILTER CELLS
# ------------------------------------------------------------------------------
md <- obj[[]]

cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    md[[cluster_col]] %in% clusters_use
]

obj <- subset(obj, cells = cells_keep)
md <- obj[[]]
expr <- expr[present_genes, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# CELL COUNTS TABLE
# ------------------------------------------------------------------------------
cell_count_table <- md %>%
  dplyr::count(.data[[sample_col]], .data[[cluster_col]], name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from = .data[[cluster_col]],
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::rename(hash.ID = 1)

write.csv(
  cell_count_table,
  file.path(outdir, "LAM_sample_by_cluster_cell_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
scale_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  z
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

mean_by_group <- function(mat, groups, ordered_levels) {
  out <- lapply(ordered_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(mat)
  colnames(out) <- ordered_levels
  out
}

# ------------------------------------------------------------------------------
# SAMPLE LOOP
# ------------------------------------------------------------------------------
samples_use <- sort(unique(md[[sample_col]]))

for (s in samples_use) {

  message("Processing sample: ", s)

  cells_s <- rownames(md)[md[[sample_col]] == s]
  if (length(cells_s) == 0) next

  md_s <- md[cells_s, , drop = FALSE]
  clusters_present <- clusters_use[clusters_use %in% unique(md_s[[cluster_col]])]

  if (length(clusters_present) == 0) next

  expr_s <- expr[, cells_s, drop = FALSE]

  avg_expr <- mean_by_group(
    mat = expr_s,
    groups = md_s[[cluster_col]],
    ordered_levels = clusters_present
  )

  zmat <- scale_rows(avg_expr)
  sample_tag <- safe_name(s)

  write.csv(
    as.data.frame(avg_expr),
    file.path(outdir, paste0(sample_tag, "_avg_expression.csv"))
  )

  write.csv(
    as.data.frame(zmat),
    file.path(outdir, paste0(sample_tag, "_row_zscore.csv"))
  )

  png(
    file.path(outdir, paste0(sample_tag, "_LAM_heatmap.png")),
    width = 2400,
    height = 1800,
    res = 300
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("LAM | ", s)
  )

  dev.off()

  pdf(
    file.path(outdir, paste0(sample_tag, "_LAM_heatmap.pdf")),
    width = 8,
    height = 5.5
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("LAM | ", s)
  )

  dev.off()
}

# ------------------------------------------------------------------------------
# RUN SUMMARY
# ------------------------------------------------------------------------------
writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present_genes)),
    paste0("Signature genes missing: ", length(missing_genes)),
    paste0("Clusters requested: ", paste(clusters_use, collapse = ", ")),
    paste0("Samples processed: ", length(samples_use))
  ),
  con = file.path(outdir, "LAM_run_summary.txt")
)

message("DONE: all sample-level LAM heatmaps saved in:")
message(outdir)

# LAM Ridge confirmation

# RIDGEPLOT: LAM signature across clusters, SAMPLE as replication unit
# Same concept as your original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/LAM_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(tibble)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "LAM_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

sig_genes <- c(
  "Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5","Lgals1","Lgals3","Cd9","Cd36"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) {
  message("Missing ignored: ", paste(missing, collapse = ", "))
}

write.csv(
  data.frame(gene = present),
  file.path(outdir, "LAM_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing),
  file.path(outdir, "LAM_missing_genes.csv"),
  row.names = FALSE
)

if (length(present) == 0) {
  stop("No LAM signature genes found in normalized RNA data.")
}

score_base <- "LAM_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

if (length(cells_keep) == 0) {
  stop("No cells left after filtering by selected clusters and valid hash.ID.")
}

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(
  plot_df,
  file.path(outdir, "LAM_sample_level_scores_by_cluster.csv"),
  row.names = FALSE
)

# optional count table
count_df <- plot_df %>%
  dplyr::count(final_clusters, name = "n_samples")

write.csv(
  count_df,
  file.path(outdir, "LAM_number_of_samples_per_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(
  plot_df,
  aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)
) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "LAM signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(
  file.path(outdir, "LAM_sample_replicate_ridgeplot.png"),
  p_ridge,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "LAM_sample_replicate_ridgeplot.pdf"),
  p_ridge,
  width = 8.5,
  height = 6.5
)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "LAM_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)

# m1PROINFLAMMAATORY

# RIDGEPLOT: LAM signature across clusters, SAMPLE as replication unit
# Same concept as your original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/LAM_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(tibble)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "LAM_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

sig_genes <- c(
  "Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5","Lgals1","Lgals3","Cd9","Cd36"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) {
  message("Missing ignored: ", paste(missing, collapse = ", "))
}

write.csv(
  data.frame(gene = present),
  file.path(outdir, "LAM_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing),
  file.path(outdir, "LAM_missing_genes.csv"),
  row.names = FALSE
)

if (length(present) == 0) {
  stop("No LAM signature genes found in normalized RNA data.")
}

score_base <- "LAM_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

if (length(cells_keep) == 0) {
  stop("No cells left after filtering by selected clusters and valid hash.ID.")
}

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(
  plot_df,
  file.path(outdir, "LAM_sample_level_scores_by_cluster.csv"),
  row.names = FALSE
)

# optional count table
count_df <- plot_df %>%
  dplyr::count(final_clusters, name = "n_samples")

write.csv(
  count_df,
  file.path(outdir, "LAM_number_of_samples_per_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(
  plot_df,
  aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)
) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "LAM signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(
  file.path(outdir, "LAM_sample_replicate_ridgeplot.png"),
  p_ridge,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "LAM_sample_replicate_ridgeplot.pdf"),
  p_ridge,
  width = 8.5,
  height = 6.5
)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "LAM_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)

# m1 proinflammatory heatmaps per sample - ccr7 replaced with cd80
# FINAL CLEAN SAMPLE-LEVEL HEATMAPS: M1 PROINFLAMMATORY SIGNATURE
# Robust version
# - forces RNA log-normalization
# - one heatmap per sample (hash.ID)
# - averages normalized expression by cluster inside each sample
# - row z-score only for visualization
# - Seurat v5 safe
# ------------------------------------------------------------------------------
# INPUT:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# OUTPUT:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/M1_PROINFLAMMATORY

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "M1_PROINFLAMMATORY"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14")

# ------------------------------------------------------------------------------
# M1 PROINFLAMMATORY SIGNATURE
# ------------------------------------------------------------------------------
sig_genes <- c(
  "Nos2","Il1a","Il1b","Tnf","Il6","Cd80","Cd86","Cxcl9","Cxcl10","Stat1"
)

# ------------------------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
stopifnot(sample_col %in% colnames(obj[[]]))

obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

# ------------------------------------------------------------------------------
# FORCE NORMALIZATION
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

present_genes <- intersect(sig_genes, rownames(expr))
missing_genes <- setdiff(sig_genes, rownames(expr))

write.csv(
  data.frame(gene = present_genes),
  file.path(outdir, "M1_PROINFLAMMATORY_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing_genes),
  file.path(outdir, "M1_PROINFLAMMATORY_missing_genes.csv"),
  row.names = FALSE
)

if (length(present_genes) == 0) {
  stop("No M1 proinflammatory genes found after normalization.")
}

# ------------------------------------------------------------------------------
# FILTER CELLS
# ------------------------------------------------------------------------------
md <- obj[[]]

cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    md[[cluster_col]] %in% clusters_use
]

obj <- subset(obj, cells = cells_keep)
md <- obj[[]]
expr <- expr[present_genes, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# CELL COUNTS TABLE
# ------------------------------------------------------------------------------
cell_count_table <- md %>%
  dplyr::count(.data[[sample_col]], .data[[cluster_col]], name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from = .data[[cluster_col]],
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::rename(hash.ID = 1)

write.csv(
  cell_count_table,
  file.path(outdir, "M1_PROINFLAMMATORY_sample_by_cluster_cell_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
scale_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  z
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

mean_by_group <- function(mat, groups, ordered_levels) {
  out <- lapply(ordered_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(mat)
  colnames(out) <- ordered_levels
  out
}

# ------------------------------------------------------------------------------
# SAMPLE LOOP
# ------------------------------------------------------------------------------
samples_use <- sort(unique(md[[sample_col]]))

for (s in samples_use) {

  message("Processing sample: ", s)

  cells_s <- rownames(md)[md[[sample_col]] == s]
  if (length(cells_s) == 0) next

  md_s <- md[cells_s, , drop = FALSE]
  clusters_present <- clusters_use[clusters_use %in% unique(md_s[[cluster_col]])]

  if (length(clusters_present) == 0) next

  expr_s <- expr[, cells_s, drop = FALSE]

  avg_expr <- mean_by_group(
    mat = expr_s,
    groups = md_s[[cluster_col]],
    ordered_levels = clusters_present
  )

  zmat <- scale_rows(avg_expr)
  sample_tag <- safe_name(s)

  write.csv(
    as.data.frame(avg_expr),
    file.path(outdir, paste0(sample_tag, "_avg_expression.csv"))
  )

  write.csv(
    as.data.frame(zmat),
    file.path(outdir, paste0(sample_tag, "_row_zscore.csv"))
  )

  png(
    file.path(outdir, paste0(sample_tag, "_M1_PROINFLAMMATORY_heatmap.png")),
    width = 2400,
    height = 1800,
    res = 300
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("M1 PROINFLAMMATORY | ", s)
  )

  dev.off()

  pdf(
    file.path(outdir, paste0(sample_tag, "_M1_PROINFLAMMATORY_heatmap.pdf")),
    width = 8,
    height = 5.5
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("M1 PROINFLAMMATORY | ", s)
  )

  dev.off()
}

# ------------------------------------------------------------------------------
# RUN SUMMARY
# ------------------------------------------------------------------------------
writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present_genes)),
    paste0("Signature genes missing: ", length(missing_genes)),
    paste0("Clusters requested: ", paste(clusters_use, collapse = ", ")),
    paste0("Samples processed: ", length(samples_use))
  ),
  con = file.path(outdir, "M1_PROINFLAMMATORY_run_summary.txt")
)

message("DONE: all sample-level M1 PROINFLAMMATORY heatmaps saved in:")
message(outdir)

# RIDGEPLOT CONFIRMATION M1 PROINFLAMMATORY ccr7 replaced for cd80
# RIDGEPLOT: M1 PROINFLAMMATORY signature across clusters, SAMPLE as replication unit
# Same concept as your original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/M1_PROINFLAMMATORY_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(tibble)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "M1_PROINFLAMMATORY_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14")

sig_genes <- c(
  "Nos2","Il1a","Il1b","Tnf","Il6","Cd80","Cd86","Cxcl9","Cxcl10","Stat1"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) {
  message("Missing ignored: ", paste(missing, collapse = ", "))
}

write.csv(
  data.frame(gene = present),
  file.path(outdir, "M1_PROINFLAMMATORY_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing),
  file.path(outdir, "M1_PROINFLAMMATORY_missing_genes.csv"),
  row.names = FALSE
)

if (length(present) == 0) {
  stop("No M1 proinflammatory signature genes found in normalized RNA data.")
}

score_base <- "M1_PROINFLAMMATORY_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

if (length(cells_keep) == 0) {
  stop("No cells left after filtering by selected clusters and valid hash.ID.")
}

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(
  plot_df,
  file.path(outdir, "M1_PROINFLAMMATORY_sample_level_scores_by_cluster.csv"),
  row.names = FALSE
)

count_df <- plot_df %>%
  dplyr::count(final_clusters, name = "n_samples")

write.csv(
  count_df,
  file.path(outdir, "M1_PROINFLAMMATORY_number_of_samples_per_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(
  plot_df,
  aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)
) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "M1 PROINFLAMMATORY signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(
  file.path(outdir, "M1_PROINFLAMMATORY_sample_replicate_ridgeplot.png"),
  p_ridge,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "M1_PROINFLAMMATORY_sample_replicate_ridgeplot.pdf"),
  p_ridge,
  width = 8.5,
  height = 6.5
)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "M1_PROINFLAMMATORY_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)

#m2 immunoregulatory HEATMAPS SAMPLE
# FINAL CLEAN SAMPLE-LEVEL HEATMAPS: M2 IMMUNOREGULATORY SIGNATURE
# Robust version
# - forces RNA log-normalization
# - one heatmap per sample (hash.ID)
# - averages normalized expression by cluster inside each sample
# - row z-score only for visualization
# - Seurat v5 safe
# ------------------------------------------------------------------------------
# INPUT:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# OUTPUT:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/M2_IMMUNOREGULATORY

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "M2_IMMUNOREGULATORY"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# LOAD
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14")

# ------------------------------------------------------------------------------
# M2 IMMUNOREGULATORY SIGNATURE
# ------------------------------------------------------------------------------
sig_genes <- c(
  "Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22","Stat6"
)

# ------------------------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
stopifnot(sample_col %in% colnames(obj[[]]))

obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

# ------------------------------------------------------------------------------
# FORCE NORMALIZATION
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

expr <- GetAssayData(obj, assay = "RNA", layer = "data")

present_genes <- intersect(sig_genes, rownames(expr))
missing_genes <- setdiff(sig_genes, rownames(expr))

write.csv(
  data.frame(gene = present_genes),
  file.path(outdir, "M2_IMMUNOREGULATORY_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing_genes),
  file.path(outdir, "M2_IMMUNOREGULATORY_missing_genes.csv"),
  row.names = FALSE
)

if (length(present_genes) == 0) {
  stop("No M2 immunoregulatory genes found after normalization.")
}

# ------------------------------------------------------------------------------
# FILTER CELLS
# ------------------------------------------------------------------------------
md <- obj[[]]

cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    md[[cluster_col]] %in% clusters_use
]

obj <- subset(obj, cells = cells_keep)
md <- obj[[]]
expr <- expr[present_genes, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# CELL COUNTS TABLE
# ------------------------------------------------------------------------------
cell_count_table <- md %>%
  dplyr::count(.data[[sample_col]], .data[[cluster_col]], name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from = .data[[cluster_col]],
    values_from = n_cells,
    values_fill = 0
  ) %>%
  dplyr::rename(hash.ID = 1)

write.csv(
  cell_count_table,
  file.path(outdir, "M2_IMMUNOREGULATORY_sample_by_cluster_cell_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
scale_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.na(z)] <- 0
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  z
}

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

mean_by_group <- function(mat, groups, ordered_levels) {
  out <- lapply(ordered_levels, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(mat)
  colnames(out) <- ordered_levels
  out
}

# ------------------------------------------------------------------------------
# SAMPLE LOOP
# ------------------------------------------------------------------------------
samples_use <- sort(unique(md[[sample_col]]))

for (s in samples_use) {

  message("Processing sample: ", s)

  cells_s <- rownames(md)[md[[sample_col]] == s]
  if (length(cells_s) == 0) next

  md_s <- md[cells_s, , drop = FALSE]
  clusters_present <- clusters_use[clusters_use %in% unique(md_s[[cluster_col]])]

  if (length(clusters_present) == 0) next

  expr_s <- expr[, cells_s, drop = FALSE]

  avg_expr <- mean_by_group(
    mat = expr_s,
    groups = md_s[[cluster_col]],
    ordered_levels = clusters_present
  )

  zmat <- scale_rows(avg_expr)
  sample_tag <- safe_name(s)

  write.csv(
    as.data.frame(avg_expr),
    file.path(outdir, paste0(sample_tag, "_avg_expression.csv"))
  )

  write.csv(
    as.data.frame(zmat),
    file.path(outdir, paste0(sample_tag, "_row_zscore.csv"))
  )

  png(
    file.path(outdir, paste0(sample_tag, "_M2_IMMUNOREGULATORY_heatmap.png")),
    width = 2400,
    height = 1800,
    res = 300
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("M2 IMMUNOREGULATORY | ", s)
  )

  dev.off()

  pdf(
    file.path(outdir, paste0(sample_tag, "_M2_IMMUNOREGULATORY_heatmap.pdf")),
    width = 8,
    height = 5.5
  )

  pheatmap(
    zmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0("M2 IMMUNOREGULATORY | ", s)
  )

  dev.off()
}

# ------------------------------------------------------------------------------
# RUN SUMMARY
# ------------------------------------------------------------------------------
writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present_genes)),
    paste0("Signature genes missing: ", length(missing_genes)),
    paste0("Clusters requested: ", paste(clusters_use, collapse = ", ")),
    paste0("Samples processed: ", length(samples_use))
  ),
  con = file.path(outdir, "M2_IMMUNOREGULATORY_run_summary.txt")
)

message("DONE: all sample-level M2 IMMUNOREGULATORY heatmaps saved in:")
message(outdir)

# RIDGEPLOT CONFIRMATION

# RIDGEPLOT: M2 IMMUNOREGULATORY signature across clusters, SAMPLE as replication unit
# Same concept as your original script, but:
#   - each point = one sample (hash.ID) within one cluster
#   - ridgeplot shows distribution of per-sample mean module scores across clusters
# ------------------------------------------------------------------------------
# Input:
#   ../Results/R_Output/seurat_final_clean4.rds
# ------------------------------------------------------------------------------
# Output:
#   output/20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION/M2_IMMUNOREGULATORY_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(tibble)
})

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
rds_path <- Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds")
output_root <- Sys.getenv("OUTPUT_DIR", unset = "output")
outdir <- file.path(
  output_root,
  "20260330 SAMPLE LEVEL SIGNATURE CONFIRMATION",
  "M2_IMMUNOREGULATORY_RIDGEPLOT_CLUSTERWISE_SAMPLE_REPLICATES"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# load
# ------------------------------------------------------------------------------
obj <- readRDS(rds_path)
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","2","3","4","5","14")

sig_genes <- c(
  "Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22","Stat6"
)

# ------------------------------------------------------------------------------
# normalize + module score
# ------------------------------------------------------------------------------
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)

message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) {
  message("Missing ignored: ", paste(missing, collapse = ", "))
}

write.csv(
  data.frame(gene = present),
  file.path(outdir, "M2_IMMUNOREGULATORY_present_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = missing),
  file.path(outdir, "M2_IMMUNOREGULATORY_missing_genes.csv"),
  row.names = FALSE
)

if (length(present) == 0) {
  stop("No M2 immunoregulatory signature genes found in normalized RNA data.")
}

score_base <- "M2_IMMUNOREGULATORY_Sig"
score_col <- paste0(score_base, "1")

obj <- AddModuleScore(
  obj,
  features = list(present),
  name = score_base,
  assay = "RNA",
  search = FALSE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# filter like original
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
obj[[sample_col]] <- as.character(obj[[sample_col]][,1])

md <- obj[[]]
cells_keep <- rownames(md)[
  !is.na(md[[sample_col]]) &
    md[[sample_col]] != "" &
    !is.na(md[[cluster_col]]) &
    md[[cluster_col]] %in% clusters_use
]

if (length(cells_keep) == 0) {
  stop("No cells left after filtering by selected clusters and valid hash.ID.")
}

obj_sub <- subset(obj, cells = cells_keep)

# ------------------------------------------------------------------------------
# collapse to sample-level replication
# ------------------------------------------------------------------------------
plot_df <- obj_sub[[]] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    hash.ID = .data[[sample_col]],
    final_clusters = as.character(.data[[cluster_col]]),
    module_score = .data[[score_col]]
  ) %>%
  dplyr::group_by(hash.ID, final_clusters) %>%
  dplyr::summarise(
    module_score = mean(module_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df$final_clusters <- factor(plot_df$final_clusters, levels = clusters_use)

write.csv(
  plot_df,
  file.path(outdir, "M2_IMMUNOREGULATORY_sample_level_scores_by_cluster.csv"),
  row.names = FALSE
)

count_df <- plot_df %>%
  dplyr::count(final_clusters, name = "n_samples")

write.csv(
  count_df,
  file.path(outdir, "M2_IMMUNOREGULATORY_number_of_samples_per_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# ridgeplot: same outcome, sample as replicate
# ------------------------------------------------------------------------------
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- ggplot(
  plot_df,
  aes(x = module_score, y = final_clusters, fill = final_clusters, color = final_clusters)
) +
  ggridges::geom_density_ridges(
    alpha = 0.55,
    linewidth = 0.25,
    scale = 1.1,
    rel_min_height = 0.001
  ) +
  scale_fill_manual(values = ridge_cols, drop = FALSE) +
  scale_color_manual(values = ridge_cols, drop = FALSE) +
  labs(
    title = "M2 IMMUNOREGULATORY signature module score across clusters",
    subtitle = "Sample-level replication unit (mean score per hash.ID within each cluster)",
    x = paste0(score_col, " (sample-level mean module score)"),
    y = "final_clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

ggsave(
  file.path(outdir, "M2_IMMUNOREGULATORY_sample_replicate_ridgeplot.png"),
  p_ridge,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "M2_IMMUNOREGULATORY_sample_replicate_ridgeplot.pdf"),
  p_ridge,
  width = 8.5,
  height = 6.5
)

writeLines(
  c(
    paste0("Input object: ", rds_path),
    paste0("Output folder: ", outdir),
    paste0("Clusters used: ", paste(clusters_use, collapse = ", ")),
    paste0("Signature genes total: ", length(sig_genes)),
    paste0("Signature genes present: ", length(present)),
    paste0("Replication unit: sample (hash.ID within cluster)"),
    paste0("Number of sample-cluster observations: ", nrow(plot_df))
  ),
  con = file.path(outdir, "M2_IMMUNOREGULATORY_sample_replicate_ridgeplot_summary.txt")
)

message("Done. Results saved in: ", outdir)
