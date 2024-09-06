rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/HmHCC/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")
load("seuList_processed.rda")

i <- as.integer(commandArgs(TRUE))
seu <- seuList[[i]]
rm(seuList)

ana_wd <- paste0(work_wd, "HCC", i, "/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

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


################################################################################
# 2. annotation
################################################################################
seu <- pdistance(seu, reduction = "fastimg", assay.name = "distce")
distce <- Seurat::GetAssayData(object = seu, slot = "data", assay = "distce")
gene.use <- rownames(distce)[apply(distce, 1, function(x) all(!is.na(x)))]

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(
    seu,
    marker.freq = marker.freq, reduction.name = "caesar", gene.use = gene.use
)

res_CelliD <- CelliD.Annotation(seu, markerList[[1]], q_est = q_est)

anno.df <- data.frame(
    res_CAESAR$CAESAR, res_CAESAR$CAESARunasg, res_CAESAR$CAESARconf,
    res_CelliD$pred, res_CelliD$pred_unassign
)
colnames(anno.df) <- c("CAESAR", "CAESARunasg", "CAESARconf", "CelliD", "CelliDunasg")
seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

ave.dist <- res_CAESAR$ave.dist
save(ave.dist, file = "ave.dist.rda")

coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 500
plot_width <- coord[1] / 500

## 2.1 plot methods - fig. 3.c
for (pred in c("CAESARunasg", "CelliDunasg")) {
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



## 2.2 measures - fig. 3.e
df <- read.csv(paste0(work_wd, "manual_labels/", "labels_HCC", i, ".csv"))
rownames(df) <- df$spot
df <- df[colnames(seu), ]

seu$manual_labels <- df$label

labels201 <- function(y) {
    y0 <- rep(1, length(y))
    y0[y == "Malignant cell"] <- 0
    y0[y == "HPC-like"] <- 0
    y0[y == "unassigned"] <- 0.5
    y0
}
y1 <- labels201(seu$CAESARunasg)
y2 <- labels201(seu$CelliDunasg)
y <- seu$manual_labels
res_acc <- c("CAESAR" = mean(y1 == y), "CELLID" = mean(y2 == y))

print(res_acc)

save(res_acc, file = "res_acc.rda")
