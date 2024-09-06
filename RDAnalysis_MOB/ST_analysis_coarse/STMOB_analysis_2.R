rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/STMOB/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")
load("seu_processed.rda")

ana_wd <- paste0(work_wd, "STMOB/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

################################################################################
# 1. co-embedding
################################################################################
## 1.1 run fast
pos <- Embeddings(seu, "pos")
seu <- CAESAR.coembedding(seu, pos, q = q_est)

## 1.2 run mca
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
seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
distce <- Seurat::GetAssayData(object = seu, slot = "data", assay = "distce")
gene.use <- rownames(distce)[apply(distce, 1, function(x) all(!is.na(x)))]

marker.freq <- markerList2mat(markerList)
res_CAESAR <- CAESAR.annotation(seu, marker.freq = marker.freq, reduction.name = "caesar")

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
plot_height <- coord[2] / 5
plot_width <- coord[1] / 5

## 2.1 plot methods - fig. 4.a
for (pred in c("CAESARunasg", "CelliDunasg")) {
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
        file = paste0(ana_wd, "DimPlot_pos_", pred, ".png"),
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
    file = paste0(ana_wd, "DimPlot_pos_manual_annotation.png"),
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
    file = paste0(ana_wd, "FeaturePlot_pos_", "CAESARcf", ".png"),
    plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
)

## 2.2 measures
acc_st <- function(manual_annotation, pred) {
    manual_annotation <- as.character(manual_annotation)
    pred <- as.character(pred)
    manual_annotation[manual_annotation == "GCL"] <- "GC"
    manual_annotation[manual_annotation == "MCL"] <- "M/TC"
    manual_annotation[manual_annotation == "ONL"] <- "OSNs"
    manual_annotation[manual_annotation == "GL"] <- "PGC"
    return(mean(manual_annotation == pred))
}

res_acc <- c(
    CAESAR = acc_st(seu$manual_annotation, seu$CAESARunasg),
    CelliD = acc_st(seu$manual_annotation, seu$CelliDunasg)
)

res_pou <- c(
    CAESAR = mean(seu$CAESARunasg == "unassigned"),
    CelliD = mean(seu$CelliDunasg == "unassigned")
)

res_asw <- c(
    CAESAR = asw(Embeddings(seu, "caesar"), seu$manual_annotation),
    CelliD = asw(Embeddings(seu, "mca"), seu$manual_annotation)
)
names(res_asw) <- c("CAESAR", "CelliD")

print(res_acc)
print(res_pou)
print(res_asw)

save(res_pou, res_acc, res_asw, file = paste0(work_wd, "res_measures.rda"))

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")


## 2.5 confusion matrix - fig. 4.b
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
    file = paste0(ana_wd, "confusion_matrix_CAESARunasg.png"),
    plot = plot, width = 5, height = 3.4, units = "in", dpi = 200
)
