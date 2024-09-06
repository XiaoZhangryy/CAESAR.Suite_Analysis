rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")
load("seu_processed.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")
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
plot_height <- coord[2] / 2000
plot_width <- coord[1] / 2000

## 2.1 plot methods - fig. 4.c right panel
for (pred in c("CAESARunasg", "CelliDunasg")) {
    Idents(seu) <- factor(seu[[pred]][, 1], levels = celltypes)
    plot <- DimPlot(seu, reduction = "pos", cols = cols, pt.size = 0.4) +
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

Idents(seu) <- factor(seu$CAESARunasg, levels = names(cols))

pred <- as.character(seu$CAESARunasg)
pred[pred == "Neuron.OSN"] <- "OSNs"
levels_pred <- names(cols)
levels_pred[levels_pred == "Neuron.OSN"] <- "OSNs"
Idents(seu) <- factor(pred, levels = levels_pred)
cols_pred <- cols
names(cols_pred) <- levels_pred

plot2 <- DimPlot(
    seu,
    reduction = "pos", cols = cols_pred, split.by = "ident",
    group.by = "ident", ncol = 3, pt.size = 0.4, raster = FALSE
) +
    theme(
        plot.title = element_blank(),
        axis.line = element_blank(), # Remove the axis lines
        axis.text.x = element_blank(), # Remove the text on the x-axis
        axis.text.y = element_blank(), # Remove the text on the y-axis
        axis.ticks = element_blank(), # Remove the ticks
        axis.title.x = element_blank(), # Remove the x-axis title
        axis.title.y = element_blank(),
        legend.position = "right",
        panel.background = element_rect(fill = "grey90", colour = "grey50")
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))
ggsave(
    file = paste0(ana_wd, "DimPlot_pos_split_CAESARunasg.png"),
    plot = plot2, width = plot_width * 2, height = plot_height * 4, units = "in", dpi = 500
)

plot <- FeaturePlot(
    seu,
    reduction = "pos", features = "CAESARconf", pt.size = 0.4,
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

## 2.2 measures
res_acc <- c(
    CAESAR = mean(seu$merged_deviAnnotate == seu$CAESARunasg),
    CelliD = mean(seu$merged_deviAnnotate == seu$CelliDunasg)
)

res_pou <- c(
    CAESAR = mean(seu$CAESARunasg == "unassigned"),
    CelliD = mean(seu$CelliDunasg == "unassigned")
)

res_asw <- c(
    CAESAR = asw(Embeddings(seu, "caesar"), seu$merged_deviAnnotate),
    CelliD = asw(Embeddings(seu, "mca"), seu$merged_deviAnnotate)
)
names(res_asw) <- c("CAESAR", "CelliD")

print(res_acc)
print(res_pou)
print(res_asw)

save(res_pou, res_acc, res_asw, file = paste0(work_wd, "res_measures.rda"))

## 2.3 sgs
Idents(seu) <- factor(seu$CAESARunasg, levels = names(cols))
sg_list <- find.sig.genes(seu, distce.assay = "distce", expr.prop.cutoff = 0.1)
save(sg_list, file = "sg_list.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")

## 2.5 plot each celltype - fig. 4.e
for (ct in c("Mes", "OEC", "Neuron.OSN", "PGC", "M/TC", "GC")) {
    # Subset the Seurat object for the current cell type
    seu_sub <- subset(seu, CAESARunasg == ct)

    # Plot the current cell type
    plot <- DimPlot(seu_sub, reduction = "pos", cols = cols[ct]) + theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank()
    )
    ggsave(
        file = paste0(ana_wd, "DimPlot_pos_", gsub("[/-]", ".", ct), ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
    )
}


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
