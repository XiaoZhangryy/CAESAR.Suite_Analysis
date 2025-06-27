rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")
load("seuList_processed.rda")
load("seu_sc.rda")

i <- as.integer(commandArgs(TRUE))
seu <- seuList[[i]]
rm(seuList)

ana_wd <- paste0(work_wd, "BC", i, "/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results/")
dir.create(annotation_wd, showWarnings = FALSE)

load_ref <- function(im) {
    load(paste0(work_wd, "seu_scList.rda"))
    seu_scList[[im]]
}

################################################################################
# 1. co-embedding
################################################################################
## 1.1 obtain image features
obatin_feature_img <- function(i) {
    img_path <- "/share/analysisdata/liuw/coembed_image/RealData/Xenium_hBreast4/processdata/"

    metafiles <- c(
        "meta_data_sample1.csv", "meta_data_sample2.csv",
        "meta_data_sample3.csv", "meta_data_sample4.csv"
    )

    featurefiles <- c(
        "crop_xenium_prerelease_jul12_hBreast_rep1.jpgfeature_img.csv",
        "crop_xenium_prerelease_jul12_hBreast_rep2.jpegfeature_img.csv",
        "crop_Rep1_Full_he_image.jpegfeature_img.csv",
        "crop_Rep2_Full_he_image.jpegfeature_img.csv"
    )

    meta <- read.csv(paste0(img_path, metafiles[i]))
    feature <- read.csv(paste0(img_path, featurefiles[i]), header = FALSE)
    rownames(feature) <- meta$X

    spots <- colnames(seu)
    if (!all(spots %in% rownames(feature))) {
        stop("img features does not match with seurat object")
    }

    feature <- feature[spots, ]
    feature
}

## 1.2 run fast_img
feature_img <- obatin_feature_img(i)
pos <- Embeddings(seu, "pos")
seu <- CAESAR.coembedding.image(seu, feature_img, pos, q = q_est)

## 1.3 run mca
seu <- RunMCA(seu, nmcs = q_est, reduction.name = "mca")


################################################################################
# 2. RCTD
################################################################################
if (file.exists(paste0(ana_wd, "RCTD_results/RCTD_results.rda"))) {
    load(paste0(ana_wd, "RCTD_results/RCTD_results.rda"))
} else {
    library(spacexr)
    library(Matrix)

    reference <- Reference(
        seu_sc@assays$RNA@counts,
        factor(seu_sc$merged_clusters, levels = celltypes),
        seu_sc$nCount_RNA
    )

    query <- SpatialRNA(
        seu@meta.data[, c("x_centroid", "y_centroid")],
        seu@assays$RNA@counts,
        seu$nCount_RNA
    )

    myRCTD <- create.RCTD(query, reference, max_cores = 8, UMI_min = 0, counts_MIN = 0)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")


    resultsdir <- paste0(ana_wd, "RCTD_results")
    dir.create(resultsdir, showWarnings = FALSE)
    setwd(resultsdir)
    save(myRCTD, file = "myRCTD.rda")

    results <- myRCTD@results
    # normalize the cell type proportions to sum to 1.
    norm_weights <- normalize_weights(results$weights)
    cell_type_names <- myRCTD@cell_type_info$info[[2]] # list of cell type names
    spatialRNA <- myRCTD@spatialRNA

    # make the plots
    # Plots the confident weights for each cell type as in full_mode (saved as
    # 'results/cell_type_weights_unthreshold.pdf')
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots all weights for each cell type as in full_mode. (saved as
    # 'results/cell_type_weights.pdf')
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots the weights for each cell type as in doublet_mode. (saved as
    # 'results/cell_type_weights_doublets.pdf')
    plot_weights_doublet(
        cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
        results$results_df
    )
    # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as
    # 'results/cell_type_occur.pdf')
    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
    # makes a map of all cell types, (saved as
    # 'results/all_cell_types.pdf')
    plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)

    RCTD_results <- results$results_df
    save(RCTD_results, file = "RCTD_results.rda")
}

seu$RCTD_first <- RCTD_results$first_type

save(seu, file = "seu.rda")

Idents(seu) <- seu$RCTD_first
degs <- FindAllMarkers(seu)
save(degs, file = "degs_RCTD_first.rda")

################################################################################
# 2. CAESAR annotation
################################################################################
setwd(ana_wd)

seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
distce <- Seurat::GetAssayData(object = seu, slot = "data", assay = "distce")

res_CelliD_List <- vector("list", length(markerList))
for (im in seq_along(markerList)) {
    marker <- list(markerList[[im]])
    sample <- names(markerList)[im]

    marker.freq <- markerList2mat(marker)
    res_CAESAR <- CAESAR.annotation(
        seu,
        marker.freq = marker.freq, reduction.name = "caesar"
    )

    anno.df <- data.frame(
        res_CAESAR$pred, res_CAESAR$pred_unassign, res_CAESAR$confidence
    )
    colnames(anno.df) <- paste0(c("CAESAR", "CAESARunasg", "CAESARconf"), "_", gsub(" ", "_", sample))

    
    result_file_name <- paste0(annotation_wd, gsub(" ", "_", sample), "_CAESAR.rda")
    save(anno.df, file = result_file_name)




    res_CelliD <- CelliD.Annotation(seu, marker[[1]])

    anno.df <- data.frame(
        res_CelliD$prediction, res_CelliD$prediction_significant, res_CelliD$significant
    )
    colnames(anno.df) <- paste0(c("CelliD", "CelliDunasg", "CelliDconf"), "_", gsub(" ", "_", sample))

    result_file_name <- paste0(annotation_wd, gsub(" ", "_", sample), "_CelliD.rda")
    save(anno.df, file = result_file_name)






    result_file_name <- paste0(annotation_wd, gsub(" ", "_", sample), "_Seurat_scmap_SingleR_scPred.rda")

    seu_sc_fov <- load_ref(im)

    
    res_Seurat <- Seurat.transfer(seu, seu_sc_fov, "merged_clusters")

    res_scmap <- scmap.transfer(seu, seu_sc_fov, "merged_clusters")

    res_SingleR <- SingleR.transfer(seu, seu_sc_fov, "merged_clusters")

    res_scPred <- scPred.transfer(seu, seu_sc_fov, "merged_clusters")


    anno.df <- data.frame(
        res_Seurat$pred, res_Seurat$pred_unassign, res_Seurat$significant,
        res_scmap$pred, res_scmap$pred_unassign, res_scmap$significant,
        res_SingleR$pred, res_SingleR$pred_unassign, res_SingleR$significant,
        res_scPred$pred, res_scPred$pred_unassign, res_scPred$significant
    )
    colnames(anno.df) <- c(
        "Seurat", "Seuratunasg", "Seuratconf",
        "scmap", "scmapunasg", "scmapconf",
        "SingleR", "SingleRunasg", "SingleRconf",
        "scPred", "scPredunasg", "scPredconf"
    )
    colnames(anno.df) <- paste0(colnames(anno.df), "_", gsub(" ", "_", sample))
    rownames(anno.df) <- colnames(seu)

    save(anno.df, file = result_file_name)
}
save(res_CelliD_List, file = "res_CelliD_List.rda")

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(seu, marker.freq = marker.freq, reduction.name = "caesar")
anno.df <- data.frame(
    res_CAESAR$pred, res_CAESAR$pred_unassign, res_CAESAR$confidence
)
colnames(anno.df) <- c("iCAESAR", "iCAESARunasg", "iCAESARconf")

result_file_name <- paste0(annotation_wd, "iCAESAR.rda")
save(anno.df, file = result_file_name)

ave.dist <- res_CAESAR$ave.dist
save(ave.dist, file = "ave.dist.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")
