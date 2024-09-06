rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")
load("seu_processed.rda")

ana_wd <- paste0(work_wd, "spATACME11/")
dir.create(ana_wd, showWarnings = FALSE)
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


################################################################################
# 1. co-embedding
################################################################################
## 1.1 run fast_img
img_path <- "/share/analysisdata/liuw/coembed_image/RealData/MouseEmbryo/ME11/"
meta <- read.csv(paste0(img_path, "meta_data.csv"))
feature_img <- fread(paste0(img_path, "feature_img.csv"))
feature_img <- as.matrix(feature_img)
rownames(feature_img) <- meta$barcode
pos <- Embeddings(seu, "pos")
seu <- CAESAR.coembedding.image(seu, feature_img, pos, q = q_est)

## 1.2 run mca
seu <- RunMCA(seu, nmcs = q_est, reduction.name = "mca")

################################################################################
# 2. annotation
################################################################################
seu <- pdistance(seu, reduction = "fastimg", assay.name = "distce")
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

## 2.1 plot methods - fig. 5.b
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
    CAESAR = asw(Embeddings(seu, "fastimg"), seu$merged_deviAnnotate),
    CelliD = asw(Embeddings(seu, "mca"), seu$merged_deviAnnotate)
)
names(res_asw) <- c("CAESAR", "CelliD")

print(res_acc)
print(res_pou)
print(res_asw)

save(res_pou, res_acc, res_asw, file = paste0(work_wd, "res_measures.rda"))

## 2.3 sgs
Idents(seu) <- factor(seu$CAESARunasg, levels = names(cols))
sg_list <- find.signature.genes.modified(seu, distce.assay = "distce", expr.prop.cutoff = 0.1)
save(sg_list, file = "sg_list.rda")

seu@assays[["distce"]] <- NULL
save(seu, file = "seu.rda")


## 2.5 plot each celltype - fig. 5.d left penal
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
