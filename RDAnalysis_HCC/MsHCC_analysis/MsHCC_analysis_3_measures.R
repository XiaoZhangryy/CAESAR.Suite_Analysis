rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/MsHCC/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")

################################################################################
# 1. measures - fig. 2.e
################################################################################
ms_res_acc <- sapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "res_acc.rda"))
    res_acc
})
colnames(ms_res_acc) <- paste0("HCC", 1:4)
save(ms_res_acc, file = "ms_res_acc.rda")
