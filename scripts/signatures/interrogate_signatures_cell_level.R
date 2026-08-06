#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/interrogate_signatures_cell_level.R
# Original file: 20260218 SIGNATURE INTERROGATION (CELL LEVEL).R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Interrogate selected immune signatures at single-cell level.
# Inputs: Canonical analysis-ready Seurat object with RNA data, UMAP, and cluster metadata.
# Outputs: Signature ridgeplots, violin plots, selected heatmaps, and feature plots.
# Assay/layer input: RNA data.
# Dependencies: Seurat, dplyr, ggplot2, ggridges, patchwork, pheatmap.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - An incomplete preserved double-ridge fragment is skipped unless prerequisite objects are available.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

ensure_rna_data <- function(object) {
  DefaultAssay(object) <- "RNA"
  rna_data <- tryCatch(
    GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(rna_data) || ncol(rna_data) == 0) {
    object <- NormalizeData(
      object,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
  object
}

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

# your signature
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
# compute module score
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)
message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "CustomSig"
score_col <- paste0(score_base, "1")

if (!score_col %in% colnames(obj[[]])) {
  obj <- AddModuleScore(obj, features = list(present), name = score_base, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# subset to requested clusters
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

Idents(obj) <- obj[[cluster_col]][,1]
obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# ridgeplot (transparent + discrete palette)
# ------------------------------------------------------------------------------
# muted, discrete, colorblind-friendly palette (11 clusters)
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- RidgePlot(
  obj_sub,
  features = score_col,
  ncol = 1
)

# adjust the ridge layer directly (Seurat returns a ggplot object)
p_ridge$layers[[1]]$aes_params$alpha <- 0.55  # more transparent
p_ridge$layers[[1]]$aes_params$size <- 0.25  # thinner outlines

p_ridge <- p_ridge +
  scale_fill_manual(values = ridge_cols) +
  scale_color_manual(values = ridge_cols) +
  labs(
    title = "Signature module score across clusters",
    subtitle = paste0("Clusters: ", paste(clusters_use, collapse = ", ")),
    x = paste0(score_col, " (module score)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_Signature_Ridgeplots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

# spranger interrogation

# RIDGEPLOT: custom signature across selected clusters (transparent + discrete palette)
# Assumes seurat_final_clean4 is already in memory
# Produces: a nice ridgeplot of the signature score across clusters:
# 0,2,3,4,5,14,13,1,7,8,9
# Styling:
#   - more transparent ridges
#   - discrete muted palette (colorblind-friendly)
#   - thin ridge outlines

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

# your signature
sig_genes <- c(
  "Ifit1bl1","Ifit3b","Ifit3","Ifit1","Cxcl10","Ifit2","Rsad2","Phf11d","Isg20","Gbp7","Slfn1","Ifi213","Fam26f","Ifi47","Ifi205","Nt5c3","Gbp2","Ifih1","Ifi206","Oasl1","Ifi211","Fgl2","Ddx58","1600014C10Rik","Slfn9","Igtp","Cmpk2","Lipg","Herc6","Ms4a6d","Ms4a4b","Endod1","Slfn8","Ifi214","Ddx60","Dhx58","Tor3a","Iigp1","A530064D06Rik","Znfx1","Axl","Ifi44","Mx1","Oas3","Gbp3","Carhsp1","Fcgr4","Gbp4","Slfn4","Treml2","Mthfr","Lpxn","Adap2","Phf11c","Svbp","AA467197","Trim30b","B430306N03Rik","Cd69","Rnpep","Ifi204","Ms4a4c","Isg15","Fcgr1","Irf7","Phf11b","Rtp4","Usp18","Ly6a","Ms4a6b","Ifi209","Mndal","Xaf1","Ccnd2","Oasl2","Zbp1","Slfn5","Ifi203","Ly6e","Samd9l","Phf11a","Ms4a6c","Slfn2","Sp100","H2-T22","Samhd1","Ctss","Ifitm3","Sat1","Ube2l6","Lgals3bp","Parp14","Pnp","Sp110","Selenow","Sdcbp","Snx2","Irgm1","Tspo","Trafd1","Stat2","Daxx","Ifi207","Fcer1g","Hck","Dck","Pttg1","Chmp4b","Rnf213","Nmi","Klrk1","Hmox2","Aftph","Lgals9","Ifi27l2a","Rnf34","Usp25","Tmem219","Wdfy1","2810474O19Rik","Ccr5","Scimp","Pstpip1","Dtx3l","Eif2ak2","Lilr4b","Ctsc","Cd86")

# ------------------------------------------------------------------------------
# compute module score
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)
message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "CustomSig"
score_col <- paste0(score_base, "1")

if (!score_col %in% colnames(obj[[]])) {
  obj <- AddModuleScore(obj, features = list(present), name = score_base, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# subset to requested clusters
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

Idents(obj) <- obj[[cluster_col]][,1]
obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# ridgeplot (transparent + discrete palette)
# ------------------------------------------------------------------------------
# muted, discrete, colorblind-friendly palette (11 clusters)
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- RidgePlot(
  obj_sub,
  features = score_col,
  ncol = 1
)

# adjust the ridge layer directly (Seurat returns a ggplot object)
p_ridge$layers[[1]]$aes_params$alpha <- 0.55  # more transparent
p_ridge$layers[[1]]$aes_params$size <- 0.25  # thinner outlines

p_ridge <- p_ridge +
  scale_fill_manual(values = ridge_cols) +
  scale_color_manual(values = ridge_cols) +
  labs(
    title = "Signature module score across clusters",
    subtitle = paste0("Clusters: ", paste(clusters_use, collapse = ", ")),
    x = paste0(score_col, " (module score)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_Signature_Ridgeplots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

# LAM SIGNATURE
# RIDGEPLOT: custom signature across selected clusters (transparent + discrete palette)
# Assumes seurat_final_clean4 is already in memory
# Produces: a nice ridgeplot of the signature score across clusters:
# 0,2,3,4,5,14,13,1,7,8,9
# Styling:
#   - more transparent ridges
#   - discrete muted palette (colorblind-friendly)
#   - thin ridge outlines

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

# your signature
sig_genes <- c("Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5","Lgals1","Lgals3","Cd9","Cd36")

# ------------------------------------------------------------------------------
# compute module score
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)
message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "CustomSig"
score_col <- paste0(score_base, "1")

if (!score_col %in% colnames(obj[[]])) {
  obj <- AddModuleScore(obj, features = list(present), name = score_base, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# subset to requested clusters
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

Idents(obj) <- obj[[cluster_col]][,1]
obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# ridgeplot (transparent + discrete palette)
# ------------------------------------------------------------------------------
# muted, discrete, colorblind-friendly palette (11 clusters)
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- RidgePlot(
  obj_sub,
  features = score_col,
  ncol = 1
)

# adjust the ridge layer directly (Seurat returns a ggplot object)
p_ridge$layers[[1]]$aes_params$alpha <- 0.55  # more transparent
p_ridge$layers[[1]]$aes_params$size <- 0.25  # thinner outlines

p_ridge <- p_ridge +
  scale_fill_manual(values = ridge_cols) +
  scale_color_manual(values = ridge_cols) +
  labs(
    title = "Signature module score across clusters",
    subtitle = paste0("Clusters: ", paste(clusters_use, collapse = ", ")),
    x = paste0(score_col, " (module score)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_Signature_Ridgeplots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

#PROINFLAMMATORY SIGNATURE

# LAM SIGNATURE
# RIDGEPLOT: custom signature across selected clusters (transparent + discrete palette)
# Assumes seurat_final_clean4 is already in memory
# Produces: a nice ridgeplot of the signature score across clusters:
# 0,2,3,4,5,14,13,1,7,8,9
# Styling:
#   - more transparent ridges
#   - discrete muted palette (colorblind-friendly)
#   - thin ridge outlines

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
})

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

# your signature
sig_genes <- c("Nos2", "Il1b", "Il1a", "Tnf", "Il6")

# ------------------------------------------------------------------------------
# compute module score
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)
message("Signature genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored (first 30): ", paste(head(missing, 30), collapse = ", "))

score_base <- "CustomSig"
score_col <- paste0(score_base, "1")

if (!score_col %in% colnames(obj[[]])) {
  obj <- AddModuleScore(obj, features = list(present), name = score_base, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# subset to requested clusters
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

Idents(obj) <- obj[[cluster_col]][,1]
obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# ridgeplot (transparent + discrete palette)
# ------------------------------------------------------------------------------
# muted, discrete, colorblind-friendly palette (11 clusters)
ridge_cols <- c(
  "#4E79A7", "#59A14F", "#9C755F", "#76B7B2", "#EDC948",
  "#B07AA1", "#FF9DA7", "#86BCB6", "#F28E2B", "#8CD17D", "#BAB0AC"
)

p_ridge <- RidgePlot(
  obj_sub,
  features = score_col,
  ncol = 1
)

# adjust the ridge layer directly (Seurat returns a ggplot object)
p_ridge$layers[[1]]$aes_params$alpha <- 0.55  # more transparent
p_ridge$layers[[1]]$aes_params$size <- 0.25  # thinner outlines

p_ridge <- p_ridge +
  scale_fill_manual(values = ridge_cols) +
  scale_color_manual(values = ridge_cols) +
  labs(
    title = "Signature module score across clusters",
    subtitle = paste0("Clusters: ", paste(clusters_use, collapse = ", ")),
    x = paste0(score_col, " (module score)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_Signature_Ridgeplots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.png"), p_ridge, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(outdir, "CustomSignature_ridgeplot.pdf"), p_ridge, width = 8.5, height = 6.5)

# INFLAMMATORY VIOLIN

# inflammatory TAM signature - minimal violin by cluster (PDF)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","13","1","7","8","9")

sig_genes <- c("Nos2","Il1a","Il1b","Tnf","Il6")  # inflammatory TAM

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_20250219_inos_plot")
outfile <- file.path(outdir, "InflammatoryTAM_signature_violin.pdf")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# module score (robust to missing genes)
# ------------------------------------------------------------------------------
present <- intersect(sig_genes, rownames(obj))
missing <- setdiff(sig_genes, present)
message("Genes present: ", length(present), "/", length(sig_genes))
if (length(missing) > 0) message("Missing ignored: ", paste(missing, collapse = ", "))

stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
Idents(obj) <- obj[[cluster_col]][,1]

obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

score_base <- "InflamTAM_sig"
score_col <- paste0(score_base, "1")
if (!score_col %in% colnames(obj_sub[[]])) {
  obj_sub <- AddModuleScore(obj_sub, features = list(present), name = score_base, verbose = FALSE)
}

# ------------------------------------------------------------------------------
# pale discrete colors (no lighten() needed)
# ------------------------------------------------------------------------------
base_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(clusters_use))
pale_cols <- sapply(base_cols, function(cl) alpha(cl, 0.55))
names(pale_cols) <- clusters_use

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------
p <- VlnPlot(
  obj_sub,
  features = score_col,
  group.by = cluster_col,
  pt.size = 0,
  cols = pale_cols,
  adjust = 1.3
) +
  geom_violin(width = 0.92, trim = TRUE, color = NA) +   # close + clean
  stat_summary(fun = median, geom = "crossbar",
               width = 0.35, fatten = 0, color = "black") +
  labs(
    title = "Inflammatory TAM signature (Nos2, Il1a, Il1b, Tnf, Il6)",
    x = "Cluster",
    y = "Module score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(p)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
pdf(outfile, width = 9, height = 4.8, onefile = TRUE)
print(p)
dev.off()
message("Saved: ", normalizePath(outfile, mustWork = FALSE))

# M1 VS M2
# M1 vs M2 (mouse) - module-score violins by cluster (minimal, pale, median)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14","1","7","8","9","13")  # REQUIRED ORDER

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_M1_VS_M2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# canonical signatures (mouse)
# ------------------------------------------------------------------------------
M1_genes <- c("Nos2","Il1a","Il1b","Tnf","Il6","Il12b","Ccr7","Cd86","Cxcl9","Cxcl10")
M2_genes <- c("Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22")

# ------------------------------------------------------------------------------
# prep identities / subset
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
Idents(obj) <- obj[[cluster_col]][,1]

obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# helper: add module score safely
# ------------------------------------------------------------------------------
add_score <- function(objx, genes, base_name) {
  present <- intersect(genes, rownames(objx))
  missing <- setdiff(genes, present)
  message(base_name, " present: ", length(present), "/", length(genes))
  if (length(missing) > 0) message(base_name, " missing ignored: ", paste(missing, collapse = ", "))
  score_col <- paste0(base_name, "1")
  if (!score_col %in% colnames(objx[[]])) {
    objx <- AddModuleScore(objx, features = list(present), name = base_name, verbose = FALSE)
  }
  list(obj = objx, score_col = score_col)
}

res1 <- add_score(obj_sub, M1_genes, "M1_sig")
obj_sub <- res1$obj; M1_col <- res1$score_col

res2 <- add_score(obj_sub, M2_genes, "M2_sig")
obj_sub <- res2$obj; M2_col <- res2$score_col

# ------------------------------------------------------------------------------
# pale per-cluster colors (no lighten() needed)
# ------------------------------------------------------------------------------
base_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(clusters_use))
pale_cols <- sapply(base_cols, function(cl) alpha(cl, 0.50))
names(pale_cols) <- clusters_use

# ------------------------------------------------------------------------------
# minimal violin function
# ------------------------------------------------------------------------------
make_vln <- function(objx, feature, title_text) {
  VlnPlot(
    objx,
    features = feature,
    group.by = cluster_col,
    pt.size = 0,
    cols = pale_cols,
    adjust = 1.3
  ) +
    geom_violin(width = 0.92, trim = TRUE, color = NA) +
    stat_summary(fun = median, geom = "crossbar",
                 width = 0.35, fatten = 0, color = "black") +
    labs(title = title_text, x = "Cluster", y = "Module score") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      legend.position = "none"
    )
}

p_M1 <- make_vln(obj_sub, M1_col, "Mouse M1 signature score")
p_M2 <- make_vln(obj_sub, M2_col, "Mouse M2 signature score")

print(p_M1)
print(p_M2)

# ------------------------------------------------------------------------------
# save (single PDF, 2 pages)
# ------------------------------------------------------------------------------
pdf(file.path(outdir, "M1_vs_M2_signature_violins.pdf"), width = 10, height = 4.8, onefile = TRUE)
print(p_M1)
print(p_M2)
dev.off()

# optional: separate PNGs
ggsave(file.path(outdir, "M1_signature_violin.png"), p_M1, width = 10, height = 4.8, dpi = 300)
ggsave(file.path(outdir, "M2_signature_violin.png"), p_M2, width = 10, height = 4.8, dpi = 300)

message("Saved to: ", normalizePath(outdir, mustWork = FALSE))

# TRIPLE VIOLIN PLOT
# TRIPLE violins (M1 vs M2 vs LAM) - less overlap + thin medians
# + dotted ref lines at max-median cluster per signature (SATURATED colors)
# + auto-separation (tiny y-nudge) if ref lines are too close
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14")

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_M1_M2_LAM")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile_pdf <- file.path(outdir, "M1_M2_LAM_violins_lessOverlap_thinMedian_refLines_separated.pdf")
outfile_png <- file.path(outdir, "M1_M2_LAM_violins_lessOverlap_thinMedian_refLines_separated.png")

# ------------------------------------------------------------------------------
# signatures
# ------------------------------------------------------------------------------
M1_genes <- c("Nos2","Il1a","Il1b","Tnf","Il6","Il12b","Ccr7","Cd86","Cxcl9","Cxcl10")
M2_genes <- c("Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22")
LAM_genes <- c("Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5","Lgals1","Lgals3","Cd9","Cd36")

# ------------------------------------------------------------------------------
# subset + order
# ------------------------------------------------------------------------------
stopifnot(cluster_col %in% colnames(obj[[]]))
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
Idents(obj) <- obj[[cluster_col]][,1]

obj_sub <- subset(obj, idents = clusters_use)
obj_sub[[cluster_col]] <- factor(obj_sub[[cluster_col]][,1], levels = clusters_use)
Idents(obj_sub) <- obj_sub[[cluster_col]][,1]

# ------------------------------------------------------------------------------
# module scores helper
# ------------------------------------------------------------------------------
add_score <- function(objx, genes, base_name) {
  present <- intersect(genes, rownames(objx))
  missing <- setdiff(genes, present)
  message(base_name, " present: ", length(present), "/", length(genes))
  if (length(missing) > 0) message(base_name, " missing ignored: ", paste(missing, collapse = ", "))

  score_col <- paste0(base_name, "1")
  if (!score_col %in% colnames(objx[[]])) {
    objx <- AddModuleScore(objx, features = list(present), name = base_name, verbose = FALSE)
  }
  list(obj = objx, score_col = score_col)
}

r1 <- add_score(obj_sub, M1_genes,  "M1_sig");  obj_sub <- r1$obj; M1_col <- r1$score_col
r2 <- add_score(obj_sub, M2_genes,  "M2_sig");  obj_sub <- r2$obj; M2_col <- r2$score_col
r3 <- add_score(obj_sub, LAM_genes, "LAM_sig"); obj_sub <- r3$obj; LAM_col <- r3$score_col

# ------------------------------------------------------------------------------
# dataframe
# ------------------------------------------------------------------------------
df <- FetchData(obj_sub, vars = c(cluster_col, M1_col, M2_col, LAM_col))
colnames(df)[colnames(df) == cluster_col] <- "cluster"

df <- df %>%
  pivot_longer(cols = c(M1_col, M2_col, LAM_col),
               names_to = "signature", values_to = "score") %>%
  mutate(
    signature = dplyr::case_when(
      signature == M1_col  ~ "M1",
      signature == M2_col  ~ "M2",
      signature == LAM_col ~ "LAM",
      TRUE ~ signature
    ),
    signature = factor(signature, levels = c("M1","M2","LAM")),
    cluster = factor(cluster, levels = clusters_use)
  )

# ------------------------------------------------------------------------------
# compute per-signature reference medians (max over clusters)
# ------------------------------------------------------------------------------
ref_lines <- df %>%
  group_by(signature, cluster) %>%
  summarise(med = median(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(signature) %>%
  slice_max(order_by = med, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(signature, best_cluster = cluster, y_true = med)

# ------------------------------------------------------------------------------
# auto-separate close lines (display-only)
# ------------------------------------------------------------------------------
# If any pair of ref lines is closer than eps, nudge them apart.
# eps is set relative to your data range so it scales nicely.
yrng <- range(df$score, na.rm = TRUE)
eps <- diff(yrng) * 0.02  # 2% of y-range (adjust if you want more/less separation)

ref_lines <- ref_lines %>%
  arrange(y_true) %>%
  mutate(y_plot = y_true)

# simple deterministic nudge for close neighbors
for (i in 2:nrow(ref_lines)) {
  if (abs(ref_lines$y_plot[i] - ref_lines$y_plot[i - 1]) < eps) {
    ref_lines$y_plot[i - 1] <- ref_lines$y_plot[i - 1] - eps/2
    ref_lines$y_plot[i] <- ref_lines$y_plot[i]     + eps/2
  }
}

print(ref_lines)  # shows best cluster + true median + plotted median

# ------------------------------------------------------------------------------
# plot controls (your "keep this" settings)
# ------------------------------------------------------------------------------
dodge <- position_dodge(width = 0.55)  # less overlap
violin_width <- 0.78
median_width <- 0.45
median_lw <- 0.55

# pale fills for violins
fill_cols_pale <- c(
  M1 = alpha("#ef4444", 0.45),
  M2 = alpha("#3b82f6", 0.45),
  LAM = alpha("#22c55e", 0.45)
)

# saturated colors for dotted ref lines
line_cols_sat <- c(
  M1 = "#dc2626",  # strong red
  M2 = "#2563eb",  # strong blue
  LAM = "#16a34a"   # strong green
)

p <- ggplot(df, aes(cluster, score, fill = signature)) +
  geom_violin(
    position = dodge,
    width = violin_width,
    trim = TRUE,
    scale = "width",
    color = "black",
    linewidth = 0.25,
    alpha = 0.50
  ) +
  # medians (thin black solid line segments)
  stat_summary(
    fun = median,
    geom = "errorbar",
    position = dodge,
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = median_width,
    linewidth = median_lw,
    color = "black"
  ) +
  # dotted ref lines (SATURATED colors; y_plot may be nudged slightly)
  geom_hline(
    data = ref_lines,
    aes(yintercept = y_plot, color = signature),
    inherit.aes = FALSE,
    linetype = 3,     # dotted, device-safe
    linewidth = 0.8
  ) +
  scale_fill_manual(values = fill_cols_pale) +
  scale_color_manual(values = line_cols_sat) +
  labs(
    title = "M1 vs M2 vs LAM module scores",
    subtitle = "Dotted lines = best-cluster medians (auto-nudged if too close)",
    x = "Cluster",
    y = "Module score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p)

# ------------------------------------------------------------------------------
# save (guaranteed)
# ------------------------------------------------------------------------------
pdf(outfile_pdf, width = 9.2, height = 4.8, onefile = TRUE)
print(p)
dev.off()
ggsave(outfile_png, plot = p, width = 9.2, height = 4.8, dpi = 300)

message("Saved: ", normalizePath(outdir, mustWork = FALSE))
print(list.files(outdir, full.names = TRUE))

# DOUBLE M1 M2 FEATURE PLOT

# ONE UMAP: M1 vs M2 signatures - ONLY selected clusters scored (others grey)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)
reduc <- "umap"

# clusters to evaluate
cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14")

outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_M1_VS_M2")
outfile <- file.path(outdir, "M1_M2_signature_state_selectedClusters_umap.pdf")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# signatures
# ------------------------------------------------------------------------------
M1_genes <- c("Nos2","Il1a","Il1b","Tnf","Il6","Il12b","Ccr7","Cd86","Cxcl9","Cxcl10")
M2_genes <- c("Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22")

thr_M1 <- 0
thr_M2 <- 0

# ------------------------------------------------------------------------------
# compute module scores
# ------------------------------------------------------------------------------
add_score <- function(objx, genes, base_name) {
  present <- intersect(genes, rownames(objx))
  score_col <- paste0(base_name, "1")
  if (!score_col %in% colnames(objx[[]])) {
    objx <- AddModuleScore(objx, features = list(present),
                           name = base_name, verbose = FALSE)
  }
  list(obj = objx, score_col = score_col)
}

r1 <- add_score(obj, M1_genes, "M1_sig"); obj <- r1$obj; M1_col <- r1$score_col
r2 <- add_score(obj, M2_genes, "M2_sig"); obj <- r2$obj; M2_col <- r2$score_col

# ------------------------------------------------------------------------------
# build state ONLY for chosen clusters
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

df <- FetchData(obj, vars = c(cluster_col, M1_col, M2_col))

state <- rep("Other", nrow(df))  # default = grey

inside <- df[[cluster_col]] %in% clusters_use

state[inside] <- dplyr::case_when(
  df[[M1_col]][inside] >  thr_M1 & df[[M2_col]][inside] <= thr_M2 ~ "M1 only",
  df[[M1_col]][inside] <= thr_M1 & df[[M2_col]][inside] >  thr_M2 ~ "M2 only",
  df[[M1_col]][inside] >  thr_M1 & df[[M2_col]][inside] >  thr_M2 ~ "Both",
  TRUE ~ "Neither"
)

obj$M1_M2_state <- factor(
  state,
  levels = c("Other","Neither","M1 only","M2 only","Both")
)

# ------------------------------------------------------------------------------
# cluster label positions
# ------------------------------------------------------------------------------
emb <- Embeddings(obj, reduction = reduc)
centers_df <- data.frame(
  x = emb[,1],
  y = emb[,2],
  cluster = as.character(Idents(obj))
) %>%
  group_by(cluster) %>%
  summarise(x = median(x), y = median(y), .groups = "drop")

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------
p <- DimPlot(
  obj,
  reduction = reduc,
  group.by = "M1_M2_state",
  cols = c(
    "Other" = "grey90",
    "Neither" = "grey70",
    "M1 only" = "#ef4444",
    "M2 only" = "#3b82f6",
    "Both" = "#8b5cf6"
  )
) +
  ggtitle("M1 vs M2 signatures (selected macrophage clusters)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 0.6),
    legend.title = element_blank()
  ) +
  geom_text(
    data = centers_df,
    aes(x = x, y = y, label = cluster),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  )

print(p)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
pdf(outfile, width = 7.5, height = 6.5)
print(p)
dev.off()

message("Saved: ", normalizePath(outfile, mustWork = FALSE))

# LAM VIOLIN AND RIDGE

# LAM signature - violin + ridgeplot
# pale multi-color palette + highest median reference line

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(scales)
  library(RColorBrewer)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14")

# ------------------------------------------------------------------------------
# output folder
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_LAM_signature_multicolor")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# LAM genes
# ------------------------------------------------------------------------------
LAM_genes <- c("Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5",
               "Lgals1","Lgals3","Cd9","Cd36")

# ------------------------------------------------------------------------------
# subset
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
Idents(obj) <- obj[[cluster_col]][,1]

obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# module score
# ------------------------------------------------------------------------------
present <- intersect(LAM_genes, rownames(obj_sub))
if (!"LAM_sig1" %in% colnames(obj_sub[[]])) {
  obj_sub <- AddModuleScore(obj_sub,
                            features = list(present),
                            name = "LAM_sig",
                            verbose = FALSE)
}
score_col <- "LAM_sig1"

# ------------------------------------------------------------------------------
# dataframe
# ------------------------------------------------------------------------------
df <- FetchData(obj_sub, vars = c(cluster_col, score_col))
colnames(df) <- c("cluster","score")
df$cluster <- factor(df$cluster, levels = clusters_use)

# ------------------------------------------------------------------------------
# highest median
# ------------------------------------------------------------------------------
highest_median <- df %>%
  group_by(cluster) %>%
  summarise(med = median(score), .groups="drop") %>%
  summarise(max(med)) %>%
  pull()

# ------------------------------------------------------------------------------
# pale multicolor palette
# ------------------------------------------------------------------------------
base_cols <- colorRampPalette(brewer.pal(8,"Set2"))(length(clusters_use))
pale_cols <- alpha(base_cols, 0.55)
names(pale_cols) <- clusters_use

# ------------------------------------------------------------------------------
# VIOLIN
# ------------------------------------------------------------------------------
p_vln <- ggplot(df, aes(cluster, score, fill = cluster)) +
  geom_violin(
    color = "black",
    linewidth = 0.25,
    width = 0.8,
    trim = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = 0.45,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_hline(
    yintercept = highest_median,
    linetype = 3,
    linewidth = 0.9,
    color = "black"
  ) +
  scale_fill_manual(values = pale_cols) +
  labs(
    title = "LAM signature score across clusters",
    subtitle = "Black dotted line = highest cluster median",
    x = "Cluster",
    y = "Module score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
  )

print(p_vln)

# ------------------------------------------------------------------------------
# RIDGEPLOT
# ------------------------------------------------------------------------------
p_ridge <- ggplot(df, aes(x = score, y = cluster, fill = cluster)) +
  geom_density_ridges(
    scale = 1.1,
    alpha = 0.7,
    color = "black",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = highest_median,
    linetype = 3,
    linewidth = 0.9,
    color = "black"
  ) +
  scale_fill_manual(values = pale_cols) +
  labs(
    title = "LAM signature distribution across clusters",
    subtitle = "Black dotted line = highest cluster median",
    x = "Module score",
    y = "Cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
pdf(file.path(outdir,"LAM_signature_violin.pdf"), width=7, height=5)
print(p_vln)
dev.off()

pdf(file.path(outdir,"LAM_signature_ridgeplot.pdf"), width=7, height=5)
print(p_ridge)
dev.off()

ggsave(file.path(outdir,"LAM_signature_violin.png"), p_vln, width=7, height=5, dpi=300)
ggsave(file.path(outdir,"LAM_signature_ridgeplot.png"), p_ridge, width=7, height=5, dpi=300)

message("Saved files:")
print(list.files(outdir, full.names = TRUE))

# OPTIMIZED RIDGEPLOT LAM

# LAM VIOLIN AND RIDGE
# LAM signature - violin + ridgeplot
# pale multi-color palette + highest median reference line

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(scales)
  library(RColorBrewer)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

cluster_col <- "final_clusters"
clusters_use <- c("0","2","3","4","5","14")

# ------------------------------------------------------------------------------
# output folder
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_LAM_signature_multicolor")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# LAM genes
# ------------------------------------------------------------------------------
LAM_genes <- c("Trem2","Lipa","Lpl","Ctsb","Ctsl","Fabp4","Fabp5",
               "Lgals1","Lgals3","Cd9","Cd36")

# ------------------------------------------------------------------------------
# subset
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
Idents(obj) <- obj[[cluster_col]][,1]

obj_sub <- subset(obj, idents = clusters_use)
Idents(obj_sub) <- factor(Idents(obj_sub), levels = clusters_use)

# ------------------------------------------------------------------------------
# module score
# ------------------------------------------------------------------------------
present <- intersect(LAM_genes, rownames(obj_sub))
if (!"LAM_sig1" %in% colnames(obj_sub[[]])) {
  obj_sub <- AddModuleScore(
    obj_sub,
    features = list(present),
    name = "LAM_sig",
    verbose = FALSE
  )
}
score_col <- "LAM_sig1"

# ------------------------------------------------------------------------------
# dataframe
# ------------------------------------------------------------------------------
df <- FetchData(obj_sub, vars = c(cluster_col, score_col))
colnames(df) <- c("cluster","score")
df$cluster <- factor(df$cluster, levels = clusters_use)

# ------------------------------------------------------------------------------
# highest median (single number)
# ------------------------------------------------------------------------------
highest_median <- df %>%
  group_by(cluster) %>%
  summarise(med = median(score), .groups = "drop") %>%
  summarise(max_med = max(med)) %>%
  pull(max_med)

# ------------------------------------------------------------------------------
# pale multicolor palette
# ------------------------------------------------------------------------------
base_cols <- colorRampPalette(brewer.pal(8,"Set2"))(length(clusters_use))
pale_cols <- alpha(base_cols, 0.40)
names(pale_cols) <- clusters_use

# ------------------------------------------------------------------------------
# VIOLIN
# ------------------------------------------------------------------------------
p_vln <- ggplot(df, aes(cluster, score, fill = cluster)) +
  geom_violin(
    color = alpha("black", 0.6),
    linewidth = 0.25,
    width = 0.85,
    trim = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = 0.45,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_hline(
    yintercept = highest_median,
    linetype = "dotted",
    linewidth = 0.9,
    color = "black"
  ) +
  scale_fill_manual(values = pale_cols) +
  labs(
    title = "LAM signature score across clusters",
    subtitle = "Dotted line = highest cluster median",
    x = "Cluster",
    y = "Module score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
  )

print(p_vln)

# ------------------------------------------------------------------------------
# RIDGEPLOT (more transparent + thin black outline + more overlap)
# ------------------------------------------------------------------------------
p_ridge <- ggplot(df, aes(x = score, y = cluster, fill = cluster)) +
  geom_density_ridges(
    scale = 2.1,              # tighter overlap (closer mountains)
    rel_min_height = 0.01,
    alpha = 0.55,             # more transparent mountains
    color = "black",          # thin circumference outline
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = highest_median,
    linetype = "dotted",
    linewidth = 0.9,
    color = "black"
  ) +
  scale_fill_manual(values = pale_cols) +
  labs(
    title = "LAM signature distribution across clusters",
    subtitle = "Vertical dotted line = highest cluster median",
    x = "Module score",
    y = "Cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    panel.spacing.y = unit(0.03, "lines") # tighter stacking
  )

print(p_ridge)

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
pdf(file.path(outdir,"LAM_signature_violin.pdf"), width=7, height=5)
print(p_vln)
dev.off()

pdf(file.path(outdir,"LAM_signature_ridgeplot.pdf"), width=7, height=5)
print(p_ridge)
dev.off()

ggsave(file.path(outdir,"LAM_signature_violin.png"), p_vln, width=7, height=5, dpi=300)
ggsave(file.path(outdir,"LAM_signature_ridgeplot.png"), p_ridge, width=7, height=5, dpi=300)

message("Saved files:")
print(list.files(outdir, full.names = TRUE))

# Double ridgeplot SIGNATURES MREG SPRANGER
# ------------------------------------------------------------------------------
# CONSTANT ALPHA + SLIGHT EMPHASIS OF TOP CLUSTER
# ------------------------------------------------------------------------------

if (!exists("df") || !("signature" %in% colnames(df)) || !exists("poly")) {
  message("Skipping double ridgeplot block: expected double-signature df/poly objects are not available in this preserved original fragment.")
} else {

sig_cols <- c("MREG signature" = "#4E79A7",
              "Spranger signature" = "#F28E2B")

base_alpha <- 0.35     # normal ridges
highlight_alpha <- 0.50  # slightly stronger for best cluster

# ------------------------------------------------------------------------------
# compute "highest median" cluster per signature
# ------------------------------------------------------------------------------
med_tbl <- df %>%
  group_by(signature, cluster) %>%
  summarise(med = median(score, na.rm = TRUE), .groups = "drop")

top_meds <- med_tbl %>%
  group_by(signature) %>%
  slice_max(order_by = med, n = 1, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------------------------
# assign alpha per ridge (cluster + signature)
# ------------------------------------------------------------------------------
alpha_tbl <- df %>%
  distinct(cluster, signature) %>%
  left_join(top_meds %>% select(signature, top_cluster = cluster),
            by = "signature") %>%
  mutate(alpha =
           ifelse(cluster == top_cluster,
                  highlight_alpha,
                  base_alpha)) %>%
  select(cluster, signature, alpha)

# attach alpha to polygon data
poly2 <- poly %>%
  left_join(alpha_tbl, by = c("cluster","signature"))

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------
p_double <- ggplot(poly2, aes(x = x)) +
  geom_polygon(
    aes(
      y = y,
      group = interaction(cluster, signature),
      fill = signature,
      color = signature,
      alpha = alpha
    ),
    linewidth = 0.25   # if error -> change to size = 0.25
  ) +
  geom_line(
    aes(
      y = y0,
      group = interaction(cluster, signature),
      color = signature,
      alpha = alpha
    ),
    linewidth = 0.25
  ) +
  geom_vline(
    data = top_meds,
    aes(xintercept = med, color = signature),
    linetype = "dotted",
    linewidth = 0.8,
    alpha = 1
  ) +
  scale_fill_manual(values = sig_cols) +
  scale_color_manual(values = sig_cols) +
  scale_alpha_identity(guide = "none") +
  scale_y_continuous(
    breaks = seq_along(clusters_use),
    labels = clusters_use,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(xlim = xlim_use) +
  labs(
    title = "Double ridgeplot: MREG vs Spranger signatures",
    subtitle = "Top cluster per signature slightly emphasized; dotted line = top median",
    x = "Module score",
    y = "Cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_double)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_20260219_double_ridge_mreg_interrogation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(outdir,
                       "double_ridge_constantAlpha_highlightTopCluster.png"),
  plot = p_double,
  width = 8.5,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = file.path(outdir,
                       "double_ridge_constantAlpha_highlightTopCluster.pdf"),
  plot = p_double,
  width = 8.5,
  height = 5.8
)

}

# DC identities feature plot

# cDC1 / cDC2 / Mature FEATURE PLOT

# ONE UMAP: cDC1 vs cDC2 vs Mature signatures - ONLY selected clusters scored (others grey)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4
obj <- ensure_rna_data(obj)

# create normalized data layer
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)

# confirm layers now include data
Layers(obj[["RNA"]])
reduc <- "umap"

# clusters to evaluate
cluster_col <- "final_clusters"
clusters_use <- c("1","7","8","9","13")

# output in NEW subfolder
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"), "13_cDC_signatures", "cDC1_cDC2_Mature_selectedClusters_allMixed")
outfile <- file.path(outdir, "cDC1_cDC2_Mature_signature_state_selectedClusters_umap.pdf")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# signatures
# ------------------------------------------------------------------------------
cDC1_genes <- c("Xcr1","Clec9a","Itgae","Cadm1", "Cd8a")
cDC2_genes <- c("Cd209a","Clec10a","Mgl2", "Sirpa")
Mature_genes <- c("Ccr7","Fscn1","Il12b","Cacnb3")

thr_cDC1 <- 0
thr_cDC2 <- 0
thr_Mature <- 0

# ------------------------------------------------------------------------------
# compute module scores
# ------------------------------------------------------------------------------
add_score <- function(objx, genes, base_name) {
  present <- intersect(genes, rownames(objx))
  score_col <- paste0(base_name, "1")
  if (!score_col %in% colnames(objx[[]])) {
    objx <- AddModuleScore(
      objx,
      features = list(present),
      name = base_name,
      verbose = FALSE
    )
  }
  list(obj = objx, score_col = score_col)
}

r1 <- add_score(obj, cDC1_genes,   "cDC1_sig");   obj <- r1$obj; cDC1_col <- r1$score_col
r2 <- add_score(obj, cDC2_genes,   "cDC2_sig");   obj <- r2$obj; cDC2_col <- r2$score_col
r3 <- add_score(obj, Mature_genes, "Mature_sig"); obj <- r3$obj; Mature_col <- r3$score_col

# ------------------------------------------------------------------------------
# build state ONLY for chosen clusters
# ------------------------------------------------------------------------------
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])
df <- FetchData(obj, vars = c(cluster_col, cDC1_col, cDC2_col, Mature_col))

state <- rep("Other", nrow(df))
inside <- df[[cluster_col]] %in% clusters_use

# boolean calls
b1 <- df[[cDC1_col]]   > thr_cDC1
b2 <- df[[cDC2_col]]   > thr_cDC2
b3 <- df[[Mature_col]] > thr_Mature

# Mixed = ANY >=2 positive programs
n_pos <- (b1 + b2 + b3)

state[inside] <- dplyr::case_when(
  n_pos[inside] == 0 ~ "Neither",
  n_pos[inside] >= 2 ~ "Mixed",
  (b1 & !b2 & !b3)[inside] ~ "cDC1 only",
  (!b1 &  b2 & !b3)[inside] ~ "cDC2 only",
  (!b1 & !b2 &  b3)[inside] ~ "Mature only",
  TRUE ~ "Neither"
)

obj$cDC_state <- factor(
  state,
  levels = c("Other","Neither","cDC1 only","cDC2 only","Mature only","Mixed")
)

# ------------------------------------------------------------------------------
# cluster label positions
# ------------------------------------------------------------------------------
emb <- Embeddings(obj, reduction = reduc)
centers_df <- data.frame(
  x = emb[,1],
  y = emb[,2],
  cluster = as.character(obj[[cluster_col]][,1])
) %>%
  group_by(cluster) %>%
  summarise(x = median(x), y = median(y), .groups = "drop")

# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------
p <- DimPlot(
  obj,
  reduction = reduc,
  group.by = "cDC_state",
  cols = c(
    "Other" = "grey90",
    "Neither" = "grey70",
    "cDC1 only" = "#ef4444",
    "cDC2 only" = "#3b82f6",
    "Mature only" = "#10b981",
    "Mixed" = "#8b5cf6"
  )
) +
  ggtitle("cDC1 vs cDC2 vs Mature signatures (clusters 1,7,8,9,13)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  geom_text(
    data = centers_df,
    aes(x = x, y = y, label = cluster),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  )

print(p)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
pdf(outfile, width = 7.5, height = 6.5)
print(p)
dev.off()

message("Saved: ", normalizePath(outfile, mustWork = FALSE))

