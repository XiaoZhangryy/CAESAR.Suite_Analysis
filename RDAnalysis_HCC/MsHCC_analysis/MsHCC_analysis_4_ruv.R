rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(Matrix)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/MsHCC/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("ms2hsgenes.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    seu
})

distList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "ave.dist.rda"))
    ave.dist
})

################################################################################
# 1. ruv
################################################################################
seuInt <- CAESAR.RUV(seuList, distList, verbose = FALSE, species = "human", custom_housekeep = NULL)

metaInt <- Reduce(rbind, lapply(seuList, function(seu) {
    as.matrix(seu@meta.data[, c("CAESAR", "CAESARunasg"), drop = FALSE])
})) %>% as.data.frame()
colnames(metaInt) <- c("cluster", "clusterua")
row.names(metaInt) <- colnames(seuInt)
seuInt <- AddMetaData(seuInt, metaInt, col.name = colnames(metaInt))
Idents(seuInt) <- factor(seuInt$clusterua, levels = names(cols))

save(seuInt, file = "seuInt.rda")





