
rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(DR.SC)
library(CoFAST)
library(gsdensity)
work_wd <- "/share/analysisdata/zhangx/CosMx_CoFAST/"
source(paste0(work_wd, "functions.R"))

sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)

i <- commandArgs(TRUE) %>% as.integer()
ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
setwd(ana_wd)
load(paste0(ana_wd, "seu_all.rda"))
load("fovs.rda")
load("q_est.rda")
load("degs.rda")

# run mca
seu <- RunMCA(seu, nmcs = q_est, reduction.name = "mca")
el <- compute.nn.edges(coembed = Embeddings(seu, "mca"))

save(seu, file = paste0(ana_wd, "seu_all.rda"))
save(el, file = paste0(ana_wd, "gsdenisty_el.rda"))

marker <- degs %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)
markerList <- split(marker$gene, marker$cluster)

## cellid
pathway_cellid <- RunCellHGT(seu, pathways = markerList, dims = 1:q_est, minSize = 2, reduction = "mca")

## gsdenisty
pathway_gsdenisty <- run.rwr.list(el = el, gene_set_list = markerList, cells = colnames(seu))

## proposed
pathway_scores <- CAESAR.enrich.score(seu, pathwaylist = markerList)

res_auc <- c(
    auc(seu$merged_clusters, pathway_scores, return.mean = TRUE),
    auc(seu$merged_clusters, t(pathway_cellid), return.mean = TRUE),
    auc(seu$merged_clusters, pathway_gsdenisty, return.mean = TRUE)
)

save(res_auc, file = paste0(ana_wd, "enrich_measures_all.rda"))
