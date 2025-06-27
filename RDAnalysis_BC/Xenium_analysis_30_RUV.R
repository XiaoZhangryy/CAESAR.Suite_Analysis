rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    seu
})

distList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "ave.dist.rda"))

    ave.dist
})

################################################################################
# 1. ruv
################################################################################
seuInt <- CAESAR.RUV(seuList, distList, verbose = FALSE, species = "human", custom_housekeep = NULL)

metaInt <- Reduce(rbind, lapply(seuList, function(seu) {
    as.matrix(seu@meta.data[, c("iCAESAR", "iCAESARunasg"), drop = FALSE])
})) %>% as.data.frame()
colnames(metaInt) <- c("cluster", "clusterua")
row.names(metaInt) <- colnames(seuInt)
seuInt <- AddMetaData(seuInt, metaInt, col.name = colnames(metaInt))
Idents(seuInt) <- factor(seuInt$clusterua, levels = names(cols))

save(seuInt, file = "seuInt.rda")
