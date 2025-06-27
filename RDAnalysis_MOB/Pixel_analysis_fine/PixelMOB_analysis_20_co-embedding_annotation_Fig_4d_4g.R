rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")

load("cols_cts_q.rda")
load("markerList.rda")
load("seu_processed.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")
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


ave.dist <- res_CAESAR$ave.dist
save(ave.dist, file = "ave.dist.rda")

coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 2000
plot_width <- coord[1] / 2000

## 2.1 plot methods
for (pred in c(
    "CAESARunasg", "CelliDunasg", "CAESAR", "CelliD",
    "Seurat", "Seuratunasg", "scmap", "scmapunasg", "SingleR", "SingleRunasg",
    "scPred", "scPredunasg"
    )) {
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
        file = paste0(figs_wd, "DimPlot_pos_", pred, ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 500
    )
}



## 2.3 sgs
Idents(seu) <- factor(seu$CAESARunasg, levels = names(cols))
seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
sg_list <- find.sig.genes(seu, distce.assay = "distce", expr.prop.cutoff = 0.1)
save(sg_list, file = "sg_list.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")


## 2.5 plot each celltype
for (ct in c("Mes", "OEC", "OSNs", "PGC", "M/TC", "GC")) {
    # Subset the Seurat object for the current cell type
    seu_sub <- subset(seu, idents = ct)

    # Plot the current cell type
    plot <- DimPlot(seu_sub, reduction = "pos", cols = cols_pred[ct]) + theme(
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
        file = paste0(figs_wd, "DimPlot_pos_", gsub("[/-]", ".", ct), ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
    )
}


################################################################################
# 3. CAESAR enrichment
################################################################################
load(paste0(work_wd, "pathwayList.rda"))

pathwaylist <- Reduce(c, pathwayList)

df_rgTest <- rgTest.pathway.seu(
    seu,
    reduction = "caesar", gene.set.list = pathwaylist,
    gene.set.cutoff = 3,
    test.type = list("ori", "gen", "wei", "max"), k = 5,
    wei.fun = c("weiMax", "weiGeo", "weiArith"),
    perm.num = 0, progress_bar = TRUE, ncores = 10, eta = 1e-4,
    genes.use = NULL, parallel = TRUE, asy_threshold = 0
)
rownames(df_rgTest) <- names(pathwaylist)
save(df_rgTest, file = "df_rgTest.rda")

pathway_scores <- pathwaylist_score(seu, pathwaylist)
save(pathway_scores, file = "pathway_scores.rda")
