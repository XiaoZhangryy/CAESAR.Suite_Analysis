rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/Xenium/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

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

################################################################################
# 2. CAESAR annotation
################################################################################
setwd(ana_wd)

CelliD.Annotation <- function(seu, markerList, q_est = 50) {
    HGTres <- RunCellHGT(seu, pathways = markerList, dims = 1:q_est, minSize = 2, reduction = "mca")
    pred <- rownames(HGTres)[apply(HGTres, 2, which.max)]
    pred_un <- apply(HGTres, 2, function(x) {
        tt <- which.max(x)
        ifelse(x[tt] > 2, yes = rownames(HGTres)[tt], no = "unassigned")
    })

    return(
        list(
            HGTres = HGTres,
            pred = pred,
            pred_unassign = pred_un
        )
    )
}

seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
distce <- Seurat::GetAssayData(object = seu, slot = "data", assay = "distce")

for (im in seq_along(markerList)) {
    marker <- list(markerList[[im]])
    sample <- names(markerList)[im]

    marker.freq <- markerList2mat(marker)
    res_CAESAR <- CAESAR.annotation(
        seu,
        marker.freq = marker.freq, reduction.name = "caesar"
    )

    res_CelliD <- CelliD.Annotation(seu, marker[[1]], q_est = q_est)

    anno.df <- data.frame(
        res_CAESAR$CAESAR, res_CAESAR$CAESARunasg, res_CAESAR$CAESARconf,
        res_CelliD$pred, res_CelliD$pred_unassign
    )
    colnames(anno.df) <- paste0(c("CAESAR", "CAESARunasg", "CAESARconf", "CelliD", "CelliDunasg"), "_", gsub(" ", "_", sample))

    seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))
}

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(seu, marker.freq = marker.freq, reduction.name = "caesar")
anno.df <- data.frame(
    res_CAESAR$CAESAR, res_CAESAR$CAESARunasg, res_CAESAR$CAESARconf
)
colnames(anno.df) <- c("iCAESAR", "iCAESARunasg", "iCAESARconf")
seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

ave.dist <- res_CAESAR$ave.dist
save(ave.dist, file = "ave.dist.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")

seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")

## 2.1 plot methods




## 2.2 measures
res_acc <- sapply(names(markerList), function(sample) {
    marker <- markerList[[sample]]

    cts_notin_sample <- setdiff(celltypes, names(marker))
    
    methods <- paste0(c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg"), "_", gsub(" ", "_", sample))
    results <- sapply(methods, function(method) {
        calculate_acc(seu@meta.data[, method], seu$RCTD_first, cts_notin_sample)
    })
    names(results) <- c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg")
    results
})

res_acc_iCAESAR <- sapply(c("iCAESAR", "iCAESARunasg"), function(method) {
    calculate_acc(seu@meta.data[, method], seu$RCTD_first, NULL)
})
names(res_acc_iCAESAR) <- c("iCAESAR", "iCAESARunasg")

save(res_acc, res_acc_iCAESAR, file = "res_acc.rda")

res_asw <- c(
    CAESAR = asw(Embeddings(seu, "caesar"), seu$RCTD_first),
    CelliD = asw(Embeddings(seu, "mca"), seu$RCTD_first)
)
save(res_asw, file = "res_asw.rda")

Idents(seu) <- seu$RCTD_first
degs <- FindAllMarkers(seu)
save(degs, file = "degs_RCTD_first.rda")

topdegs <- degs %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

res_SigScore <- c(
    CAESAR = coembed_score(seu, "caesar", "RCTD_first", topdegs),
    CelliD = coembed_score(seu, "mca", "RCTD_first", topdegs)
)
save(res_SigScore, file = "res_SigScore.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")

################################################################################
# 3. CAESAR enrichment
################################################################################
load(paste0(work_wd, "pathwayList.rda"))

pathwaylist <- Reduce(c, pathwayList)

df_rgTest <- CAESAR.enrich.pathway(
    seu, pathway.list = pathwaylist, reduction = "caesar"
)
rownames(df_rgTest) <- names(pathwaylist)
save(df_rgTest, file = "df_rgTest.rda")

pathway_scores <- CAESAR.enrich.score(seu, pathwaylist)
save(pathway_scores, file = "pathway_scores.rda")
