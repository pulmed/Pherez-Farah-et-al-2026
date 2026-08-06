#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: scripts/signatures/compare_m1_m2_signatures_with_violins.R
# Original file: 20260225 M1 VS M2 VIOLINS.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Compare M1 and M2 signatures with violin plots and statistical summaries.
# Inputs: Canonical analysis-ready Seurat object with expression data and cluster metadata.
# Outputs: M1/M2 violin PDF and statistics workbook.
# Assay/layer input: Configured expression assay from the canonical object.
# Dependencies: Seurat, dplyr, ggplot2, openxlsx.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Preserves legacy M1_M2_cluster0_vs5 output names.
# ------------------------------------------------------------------------------

if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

if (!exists("seurat_final_clean4")) {
  seurat_final_clean4 <- readRDS(Sys.getenv("SEURAT_RDS", unset = "data/seurat.rds"))
}
obj <- seurat_final_clean4

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
assay_use <- "RNA"
cluster_col <- "final_clusters"
sample_col <- "hash.ID"
clusters_use <- c("0","5")

# NEW DIRECTORY
outdir <- file.path(Sys.getenv("OUTPUT_DIR", unset = "output"),
                    "14_M1_VS_M2",
                    "Cluster0_vs5_FINAL_withStars")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pdf_out <- file.path(outdir,"M1_M2_cluster0_vs5_violin.pdf")
xlsx_out <- file.path(outdir,"Stats_cluster0_vs5.xlsx")

# signatures
M1_genes <- c("Nos2","Il1a","Il1b","Tnf","Il6","Il12b","Cd80","Cd86","Cxcl9","Cxcl10")
M2_genes <- c("Arg1","Mrc1","Retnla","Chil3","Il10","Tgfb1","Cd163","Ccl17","Ccl22")

# colors
cols_fill <- c(
  "M1_0"="#ef4444",
  "M1_5"="#f87171",
  "M2_0"="#3b82f6",
  "M2_5"="#93c5fd"
)

median_half_len <- 0.12

DefaultAssay(obj) <- assay_use
obj[[cluster_col]] <- as.character(obj[[cluster_col]][,1])

# ------------------------------------------------------------------------------
# NORMALIZATION CHECK
# ------------------------------------------------------------------------------
if (!"data" %in% SeuratObject::Layers(obj[[assay_use]])) {
  obj <- NormalizeData(obj, assay=assay_use, verbose=FALSE)
}

# ------------------------------------------------------------------------------
# MODULE SCORES
# ------------------------------------------------------------------------------
add_score <- function(objx, genes, name){
  present <- intersect(genes, rownames(objx))
  objx <- AddModuleScore(objx,
                         features=list(present),
                         name=name,
                         layer="data",
                         verbose=FALSE)
  objx
}

obj <- add_score(obj, M1_genes,"M1_sig")
obj <- add_score(obj, M2_genes,"M2_sig")

M1_col <- "M1_sig1"
M2_col <- "M2_sig1"

# ------------------------------------------------------------------------------
# BUILD LONG DATA
# ------------------------------------------------------------------------------
df <- FetchData(obj,
                vars=c(cluster_col,sample_col,M1_col,M2_col)) %>%
  filter(.data[[cluster_col]] %in% clusters_use) %>%
  mutate(
    Cluster=factor(.data[[cluster_col]],levels=clusters_use),
    hash.ID=.data[[sample_col]]
  ) %>%
  pivot_longer(cols=c(M1_col,M2_col),
               names_to="Signature",
               values_to="Score") %>%
  mutate(
    Signature=recode(Signature,
                     M1_sig1="M1",
                     M2_sig1="M2"),
    SigClust=paste0(Signature,"_",Cluster)
  )

# ------------------------------------------------------------------------------
# CELL LEVEL STATS
# ------------------------------------------------------------------------------
cell_stats <- df %>%
  group_by(Signature) %>%
  summarise(
    p_value = wilcox.test(Score~Cluster)$p.value,
    .groups="drop"
  ) %>%
  mutate(
    stars = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )

# ------------------------------------------------------------------------------
# SAMPLE LEVEL (UNPAIRED hash.ID)
# ------------------------------------------------------------------------------
per_hash <- df %>%
  group_by(hash.ID,Signature,Cluster) %>%
  summarise(mean_score=mean(Score), .groups="drop")

hash_stats <- per_hash %>%
  group_by(Signature) %>%
  summarise(
    p_value = wilcox.test(mean_score~Cluster)$p.value,
    .groups="drop"
  )

# ------------------------------------------------------------------------------
# MEDIAN LINES
# ------------------------------------------------------------------------------
med_df <- df %>%
  group_by(Signature,Cluster) %>%
  summarise(med=median(Score), .groups="drop") %>%
  mutate(
    xnum=as.numeric(Cluster),
    xstart=xnum-median_half_len,
    xend=xnum+median_half_len,
    SigClust=paste0(Signature,"_",Cluster)
  )

# star position
star_df <- df %>%
  group_by(Signature) %>%
  summarise(y=max(Score)*1.05,.groups="drop") %>%
  left_join(cell_stats,by="Signature") %>%
  mutate(x=1.5)

# ------------------------------------------------------------------------------
# PLOT
# ------------------------------------------------------------------------------
p <- ggplot(df,
            aes(x=Cluster,y=Score,
                fill=SigClust,color=SigClust))+

  geom_violin(width=.95,alpha=.65,linewidth=.6,trim=TRUE)+

  geom_segment(data=med_df,
               aes(x=xstart,xend=xend,y=med,yend=med,color=SigClust),
               inherit.aes=FALSE,
               linetype="dotted",
               linewidth=1.1)+

  geom_text(data=star_df,
            aes(x=x,y=y,label=stars),
            inherit.aes=FALSE,
            size=6,
            fontface="bold")+

  facet_wrap(~Signature,nrow=1,scales="free_y")+

  scale_fill_manual(values=cols_fill)+
  scale_color_manual(values=cols_fill)+

  labs(title="M1 and M2 module scores: Cluster 0 vs 5",
       x="Cluster",y="Module score")+

  theme_classic(base_size=12)+
  theme(
    legend.position="none",
    strip.text=element_text(face="bold"),
    panel.border=element_rect(fill=NA,linewidth=.6),
    plot.title=element_text(face="bold",hjust=.5)
  )

pdf(pdf_out,width=7.2,height=4.5)
print(p)
dev.off()

# ------------------------------------------------------------------------------
# EXPORT EXCEL
# ------------------------------------------------------------------------------
if (!requireNamespace("writexl",quietly=TRUE))
  install.packages("writexl")

writexl::write_xlsx(
  list(
    cell_level_stats=cell_stats,
    hashID_level_stats=hash_stats
  ),
  path=xlsx_out
)

message("DONE - outputs saved in:")
message(outdir)

