rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMx/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
setwd(work_wd)

assig <- function(number) {
    ns <- c(30, 29, 28, 20)
    i <- 1
    while (length(ns) > 0) {
        if (number > ns[1]) {
            number <- number - ns[1]
            ns <- ns[-1]
            i <- i + 1
        } else {
            load(paste0("/share/analysisdata/zhangx/CosMx_CoFAST/processed_data_CosMx", i, "/fovs.rda"))
            j <- as.integer(fovs[number])
            return(c(i, j))
        }
    }
}

number <- commandArgs(TRUE) %>%
    as.integer() %>%
    assig()
i <- number[1]
fov <- number[2]

ana_folder <- paste0("processed_data_CosMx", i)
ana_wd <- paste0(work_wd, ana_folder, "/")
setwd(ana_wd)

seu_name <- paste0(ana_wd, "seu_fov", fov, ".rda")
load(seu_name)

## reference data
ref_sample <- "Lung5_Rep1"
ref_i <- 1
load(paste0("/share/analysisdata/zhangx/CosMx_CoFAST/processed_data_CosMx", ref_i, "/fovs.rda"))
ref_process_wd <- paste0(work_wd, "processed_data_CosMx", ref_i, "/")
load(paste0(ref_process_wd, "q_est.rda"))
load(paste0(ref_process_wd, "cols.rda"))
cols0 <- cols

# remove genes with zero variance
nz_gene <- names(which(apply(seu@assays$RNA@counts, 1, function(x) {
    any(x > 0)
})))
seu <- seu[nz_gene, ]
rm(nz_gene)

################################################################################
# 1. annotate
################################################################################
for (rfov in fovs) {
    load(paste0(ref_process_wd, "sig_list_scRNA_fov", rfov, "_bin5.rda"))
    markerList <- lapply(sig_list, function(sig) {
        marker <- marker.select(sig, overlap.max = 2)
        overlap.max <- 2
        while (length(marker[[1]]) <= 1) {
            overlap.max <- overlap.max + 1
            marker <- marker.select(sig, overlap.max = overlap.max)
        }
        marker
    })
    n.marker <- sapply(markerList, function(marker) length(marker[[1]]))
    message("The number of markers for each reference dataset are (", paste0(n.marker, collapse = ", "), ")")

    marker.freq <- markerList2mat(markerList)

    res_CAESAR <- CAESAR.annotation(
        seu,
        marker.freq = marker.freq, reduction.name = "caesar"
    )
    res_CelliD <- CelliD.Annotation(seu, marker[[1]], q_est = q_est)

    anno.df <- data.frame(
        res_CAESAR$CAESAR, res_CAESAR$CAESARunasg, res_CAESAR$CAESARconf,
        res_CelliD$pred, res_CelliD$pred_unassign
    )
    colnames(anno.df) <- c("CAESAR", "CAESARunasg", "CAESARconf", "CelliD", "CelliDunasg")
    colnames(anno.df) <- paste0(colnames(anno.df), "_ref", ref_i, "fov", rfov)

    seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))
}

sig_list <- lapply(fovs, function(rfov) {
    load(paste0(ref_process_wd, "sig_list_scRNA_fov", rfov, "_bin5.rda"))
    sig_list[[1]]
})
markerList <- lapply(sig_list, function(sig) {
    marker <- marker.select(sig, overlap.max = 2)
    overlap.max <- 2
    while (length(marker[[1]]) <= 1) {
        overlap.max <- overlap.max + 1
        marker <- marker.select(sig, overlap.max = overlap.max)
    }
    marker
})

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(
    seu,
    marker.freq = marker.freq, reduction.name = "caesar"
)
anno.df <- data.frame(
    res_CAESAR$CAESAR, res_CAESAR$CAESARunasg, res_CAESAR$CAESARconf
)
colnames(anno.df) <- c("iCAESAR", "iCAESARunasg", "iCAESARconf")
colnames(anno.df) <- paste0(colnames(anno.df), "_ref", ref_i)

seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

################################################################################
# 2. measures
################################################################################
mymeasures <- function(y1, y2) {
    y1 <- as.character(y1)
    y2 <- as.character(y2)

    # ratio of unassigned
    ROU <- mean(y1 == "unassigned")

    # acc
    ACC <- mean(y1 == y2)

    return(setNames(c(ROU, ACC), c("ROU", "ACC")))
}

measures <- lapply(fovs, function(rfov) {
    methods <- paste0(
        c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg"), "_ref", ref_i, "fov", rfov
    )
    results <- sapply(methods, function(method) {
        mymeasures(seu@meta.data[, method], seu$merged_clusters)
    })
    colnames(results) <- c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg")
    results
})

methods <- paste0(
    c("iCAESAR", "iCAESARunasg"), "_ref", ref_i
)
iCAESAR_measures <- sapply(methods, function(method) {
    mymeasures(seu@meta.data[, method], seu$merged_clusters)
})

save(measures, iCAESAR_measures, file = paste0(ana_wd, "measures_fov", fov, "_rbin5.rda"))

# asw and sigscore are not calculated, as this is same as scenario 1







