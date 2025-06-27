rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/HmHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

load("cols_cts_q.rda")
load("markerList.rda")
load("seuList_processed.rda")

i <- as.integer(commandArgs(TRUE))
seu <- seuList[[i]]
rm(seuList)

ana_wd <- paste0(work_wd, "HCC", i, "/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)


annotation_wd <- paste0(ana_wd, "annotation_results/")
dir.create(annotation_wd, showWarnings = FALSE)

################################################################################
# 1. co-embedding
################################################################################
## 1.1 obtain image features
obatin_feature_img <- function(i) {
    img_path <- paste0("/share/analysisdata/liuw/coembed_image/RealData/HCC4/HCC", i, "/")

    meta <- read.csv(paste0(img_path, "subset_meta.csv"))
    feature <- read.csv(paste0(img_path, "feature_img.csv"), header = FALSE)
    rownames(feature) <- meta$barcode

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
# 2. annotation
################################################################################
seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
distce <- Seurat::GetAssayData(object = seu, slot = "data", assay = "distce")
gene.use <- rownames(distce)[apply(distce, 1, function(x) all(!is.na(x)))]

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(
    seu,
    marker.freq = marker.freq, reduction.name = "caesar", gene.use = gene.use
)

res_CelliD <- CelliD.Annotation(seu, markerList[[1]])

res_Seurat <- Seurat.transfer(seu, seu_sc, "merged_clusters")

res_scmap <- scmap.transfer(seu, seu_sc, "merged_clusters")

res_SingleR <- SingleR.transfer(seu, seu_sc, "merged_clusters")

res_scPred <- scPred.transfer(seu, seu_sc, "merged_clusters")


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
rownames(anno.df) <- colnames(seu)

result_file_name <- paste0(annotation_wd, "CAESAR_Seurat_scmap_SingleR_scPred_CelliD.rda")
save(anno.df, file = result_file_name)


seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

ave.dist <- res_CAESAR$ave.dist
save(ave.dist, file = "ave.dist.rda")

coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 500
plot_width <- coord[1] / 500

## 2.1 plot methods - fig. 3.c
for (pred in c(
    "CAESAR", "CAESARunasg", "CelliD", "CelliDunasg",
    "Seurat", "Seuratunasg", "scmap", "scmapunasg", "SingleR", "SingleRunasg",
    "scPred", "scPredunasg"
    )) {
    Idents(seu) <- factor(seu[[pred]][, 1], levels = celltypes)
    plot <- DimPlot(seu, reduction = "pos", cols = cols, pt.size = 1.9) +
        theme(
            axis.line = element_blank(), # Remove the axis lines
            axis.text.x = element_blank(), # Remove the text on the x-axis
            axis.text.y = element_blank(), # Remove the text on the y-axis
            axis.ticks = element_blank(), # Remove the ticks
            axis.title.x = element_blank(), # Remove the x-axis title
            axis.title.y = element_blank(),
            legend.position = "none"
        )
    ggsave(
        file = paste0(ana_wd, "DimPlot_pos_", pred, ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
    )
}

plot <- FeaturePlot(
    seu,
    reduction = "pos", features = "CAESARconf", pt.size = 1.9,
    cols = c("blue", "lightgrey"), min.cutoff = 0.0, max.cutoff = 1.0
) + theme(
    axis.line = element_blank(), # Remove the axis lines
    axis.text.x = element_blank(), # Remove the text on the x-axis
    axis.text.y = element_blank(), # Remove the text on the y-axis
    axis.ticks = element_blank(), # Remove the ticks
    axis.title.x = element_blank(), # Remove the x-axis title
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
)
ggsave(
    file = paste0(ana_wd, "FeaturePlot_pos_", "CAESARconf", ".png"),
    plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
)


df <- read.csv(paste0(work_wd, "manual_labels/", "labels_HCC", i, ".csv"))
rownames(df) <- df$spot
df <- df[colnames(seu), ]

seu$manual_labels <- df$label

save(seu, file = "seu.rda")

## 2.2 sgs
Idents(seu) <- factor(seu$CAESARunasg, levels = names(cols))
sg_list <- find.signature.genes.modified(seu, distce.assay = "distce", expr.prop.cutoff = 0.1)
save(sg_list, file = "sg_list.rda")
