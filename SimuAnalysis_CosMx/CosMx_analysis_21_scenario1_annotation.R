rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMxSimu/"
setwd(work_wd)
source("functions.R")

number <- commandArgs(TRUE) %>%
    as.integer() %>%
    assign4CosMx()
i <- number[1]
fov <- number[2]

figs_wd <- paste0(work_wd, "Figs/")
ana_folder <- paste0("processed_data_CosMx", i)
ana_wd <- paste0(work_wd, ana_folder, "/fov", fov, "/")
setwd(ana_wd)

seu_name <- paste0(ana_wd, "seu_fov", fov, ".rda")
load(seu_name)

seu_name_scenario1 <- paste0(ana_wd, "seu_fov", fov, "_scenario1.rda")

## reference data
ref_sample <- "Lung5_Rep1"
ref_i <- 1
ref_process_wd <- paste0(work_wd, "processed_data_CosMx", ref_i, "/")
load(paste0(ref_process_wd, "fovs.rda"))
load(paste0(ref_process_wd, "q_est.rda"))
load(paste0(ref_process_wd, "cols.rda"))
cols0 <- cols

# remove genes with zero variance
nz_gene <- names(which(apply(seu@assays$RNA@counts, 1, function(x) {
    any(x > 0)
})))
seu <- seu[nz_gene, ]
rm(nz_gene)

load_ref <- function(seu_name) {
    load(seu_name)
    seu
}

################################################################################
# 1. annotate
################################################################################
markerList_all <- list()
cts_refList <- setNames(vector("list", length(fovs)), fovs)
for (rfov in fovs) {
    load(paste0(ref_process_wd, "fov", rfov, "/sig_list_scRNA_fov", rfov, ".rda"))
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

    markerList_all <- c(markerList_all, markerList)
    cts_refList[[rfov]] <- names(markerList[[1]])

    marker.freq <- markerList2mat(markerList)

    res_CAESAR <- CAESAR.annotation(
        seu,
        marker.freq = marker.freq, reduction.name = "caesar"
    )
    res_CelliD <- CelliD.Annotation(seu, markerList[[1]])

    seu_sc_fov <- load_ref(
        paste0(ref_process_wd, "fov", rfov, "/seu_fov", rfov, ".rda")
    )

    res_Seurat <- Seurat.transfer(seu, seu_sc_fov, "merged_clusters")

    res_scmap <- scmap.transfer(seu, seu_sc_fov, "merged_clusters")

    res_SingleR <- SingleR.transfer(seu, seu_sc_fov, "merged_clusters")

    res_scPred <- scPred.transfer(seu, seu_sc_fov, "merged_clusters")



    result_file_name <- paste0(
        annotation_wd, "ref", ref_i, "fov", rfov,
        "_CAESAR_Seurat_scmap_SingleR_scPred_CelliD.rda")

    anno.df <- data.frame(
        res_CAESAR$pred, res_CAESAR$pred_unassign, res_CAESAR$confidence,
        res_Seurat$pred, res_Seurat$pred_unassign, res_Seurat$significant,
        res_scmap$pred, res_scmap$pred_unassign, res_scmap$significant,
        res_SingleR$pred, res_SingleR$pred_unassign, res_SingleR$significant,
        res_scPred$pred, res_scPred$pred_unassign, res_scPred$significant,
        res_CelliD$pred, res_CelliD$pred_unassign, res_CelliD$significant
    )
    colnames(anno.df) <- c(
        "CAESAR", "CAESARunasg", "CAESARconf",
        "Seurat", "Seuratunasg", "Seuratconf",
        "scmap", "scmapunasg", "scmapconf",
        "SingleR", "SingleRunasg", "SingleRconf",
        "scPred", "scPredunasg", "scPredconf",
        "CelliD", "CelliDunasg", "CelliDconf"
    )
    colnames(anno.df) <- paste0(colnames(anno.df), "_ref", ref_i, "fov", rfov)
    rownames(anno.df) <- colnames(seu)

    save(anno.df, file = result_file_name)

}

marker.freq <- markerList2mat(markerList_all)
res_CAESAR <- CAESAR.annotation(
    seu,
    marker.freq = marker.freq, reduction.name = "caesar"
)
anno.df <- data.frame(
    res_CAESAR$pred, res_CAESAR$pred_unassign, res_CAESAR$confidence
)
colnames(anno.df) <- c("iCAESAR", "iCAESARunasg", "iCAESARconf")
colnames(anno.df) <- paste0(colnames(anno.df), "_ref", ref_i)

seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

save(seu, file = seu_name_scenario1)
