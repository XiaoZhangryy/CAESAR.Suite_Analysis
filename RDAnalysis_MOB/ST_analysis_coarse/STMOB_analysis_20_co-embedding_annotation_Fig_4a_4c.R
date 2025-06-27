rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/STMOB_coarse/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")
load("markerList.rda")
load("seu_processed.rda")

ana_wd <- paste0(work_wd, "STMOB/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results/")
dir.create(annotation_wd, showWarnings = FALSE)

################################################################################
# 1. co-embedding
################################################################################
## 1.1 run caesar
pos <- Embeddings(seu, "pos")
seu <- CAESAR.coembedding(seu, pos, q = q_est)

## 1.2 run mca
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
plot_height <- coord[2] / 5
plot_width <- coord[1] / 5

## 2.1 plot methods
for (pred in c(
    "CAESARunasg", "CelliDunasg", "CAESAR", "CelliD",
    "Seurat", "Seuratunasg", "scmap", "scmapunasg", "SingleR", "SingleRunasg",
    "scPred", "scPredunasg"
    )) {
    Idents(seu) <- factor(seu[[pred]][, 1], levels = celltypes)
    plot <- DimPlot(seu, reduction = "pos", cols = cols, pt.size = 4.2) +
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
        file = paste0(figs_wd, "DimPlot_pos_", pred, ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
    )
}

Idents(seu) <- factor(seu$manual_annotation, levels = names(cols_manual))
plot <- DimPlot(seu, reduction = "pos", cols = cols_manual, pt.size = 4.2) +
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
    file = paste0(figs_wd, "DimPlot_pos_manual_annotation.png"),
    plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
)


plot <- FeaturePlot(
    seu,
    reduction = "pos", features = "CAESARconf", pt.size = 4.2,
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
    file = paste0(figs_wd, "FeaturePlot_pos_", "CAESARconf", ".png"),
    plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
)



seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")



## 2.5 confusion matrix
library(reshape2)

data <- table(seu$manual_annotation, seu$CAESARunasg)[
    c("GCL", "ONL", "GL", "MCL"), c("GC", "OSNs", "PGC", "M/TC")
]
data <- t(apply(data, 1, function(x) {
    x / sum(x)
}))
data_melt <- melt(data)
colnames(data_melt) <- c("Row", "Column", "Value")
data_melt$Row <- factor(data_melt$Row, levels = rev(rownames(data)))
data_melt$Column <- factor(data_melt$Column, levels = colnames(data))

plot <- ggplot(data_melt, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
        colors = c("#ffffff", "#ffffc0", "#ffcc00", "#ff8800", "#ff4400", "#ff0000"),
        values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
        name = "Percentage"
    ) +
    theme_minimal() +
    theme(
        # axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10), size = 18),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10), size = 18),
        plot.title = element_blank(),
        legend.title = element_text(size = 12, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold")
    ) +
    labs(x = "CAESAR predictions", y = "Manual annotations")
ggsave(
    file = paste0(figs_wd, "confusion_matrix_CAESARunasg.png"),
    plot = plot, width = 5, height = 3.4, units = "in", dpi = 200
)
