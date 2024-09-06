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

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    seu
})

################################################################################
# 1. find int sg
################################################################################
sg_List <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    if (file.exists(paste0(ana_wd, "sg_list.rda"))) {
        load(paste0(ana_wd, "sg_list.rda"))
    } else {
        seu <- seuList[[i]]
        Idents(seu) <- factor(seu$iCAESARunasg, levels = names(cols))
        seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
        sg_list <- find.sig.genes(seu, distce.assay = "distce", expr.prop.cutoff = 0.1)
        save(sg_list, file = paste0(ana_wd, "sg_list.rda"))
    }
    sg_list
})

save(sg_List, file = "sg_List.rda")

ct_ratio <- lapply(seuList, function(seu) {
    counts <- table(seu$iCAESARunasg)
    counts / sum(counts)
})

save(ct_ratio, file = "ct_ratio.rda")

sgInt <- Intsg(sg_List, 50, ct_ratio, ratio_lower_bound = 5e-3)
save(sgInt, file = "sgInt.rda")

## 1.1 plot top5 sg ont seuInt
load("seuInt.rda")

celltypes_plot_int <- celltypes[c(1, 9, 7, 5, 2, 4, 8, 6, 3)]
cols_plot_int <- cols[c(1, 9, 7, 5, 2, 4, 8, 6, 3, 10)]

library(tidyverse)
sgInt <- lapply(sgInt, head, n = 5)

sg_features <- sgInt[celltypes_plot_int] %>%
    enframe(name = "cluster", value = "gene") %>%
    unnest(gene) %>%
    .$gene %>%
    unique()
sg_features[sg_features == "C5orf46"] <- "C5ORF46"

Idents(seuInt) <- factor(seuInt$clusterua, levels = names(cols_plot_int))
newcluster.ids <- setNames(
    c(
        "B-cells", "CAFs", "Cancer epithelial", "Endothelial", "Myeloid",
        "Normal epithelial", "Plasmablasts", "PVL", "T-cells", "unassigned"
    ),
    c(
        "B-cells", "CAFs", "Cancer Epithelial", "Endothelial", "Myeloid",
        "Normal Epithelial", "Plasmablasts", "PVL", "T-cells", "unassigned"
    )
)[celltypes_plot_int]
seuInt <- RenameIdents(seuInt, newcluster.ids)
plot <- DotPlot(
    seuInt,
    # idents = celltypes_plot_int,
    idents = newcluster.ids,
    col.min = -1, col.max = 2, dot.scale = 7,
    features = sg_features, scale.min = 20, scale.max = 100
) + # NoLegend() +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 60, vjust = 1, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        # legend.box = "horizontal",
        legend.position = "right",
        legend.box.margin = margin(t = 2.5, unit = "cm"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        # legend.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.spacing.y = unit(1.5, "cm")
    )
ggsave(
    file = paste0(figs_wd, "DotPlot_top5_sigInt_CAESARunasg.png"),
    plot = plot, width = 14, height = 5, units = "in", dpi = 200
)
