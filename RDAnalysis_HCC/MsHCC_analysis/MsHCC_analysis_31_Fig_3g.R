rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/MsHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")
load("ms2hsgenes.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    seu
})

################################################################################
# 1. find sg
################################################################################
sg_List <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "sg_list.rda"))
    sg_list
})

save(sg_List, file = "sg_List.rda")

ct_ratio <- lapply(seuList, function(seu) {
    counts <- table(seu$CAESARunasg)
    counts / sum(counts)
})

save(ct_ratio, file = "ct_ratio.rda")

sgInt <- Intsg(sg_List, 50, ct_ratio, ratio_lower_bound = 5e-3)
save(sgInt, file = "sgInt.rda")

## 1.1 plot top5 sg ont seuInt
load("seuInt.rda")
load("sgInt.rda")

celltypes_plot_int <- celltypes[c(1, 5, 6, 2, 7, 3, 4)]
cols_plot_int <- cols[c(1, 5, 6, 2, 7, 3, 4, 8)]

library(tidyverse)
sgInt <- lapply(sgInt, head, n = 6)

sg_features <- sgInt[celltypes_plot_int] %>%
    enframe(name = "cluster", value = "gene") %>%
    unnest(gene) %>%
    .$gene %>%
    unique()

Idents(seuInt) <- factor(seuInt$clusterua, levels = names(cols_plot_int))

plot <- DotPlot(
    seuInt,
    idents = celltypes_plot_int, col.min = -1, col.max = 1.5, dot.scale = 7,
    features = sg_features, scale.min = 20, scale.max = 100
) + # NoLegend() +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 1, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        # legend.box = "horizontal",
        legend.position = "right",
        legend.box.margin = margin(t = 2.5, unit = "cm"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"),
        # legend.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.spacing.y = unit(-8, "cm"),
        legend.margin = margin(t = 10, unit = "cm")
    ) +
    coord_flip()
ggsave(
    file = paste0(figs_wd, "DotPlot_top5_sigInt_CAESARunasg.png"),
    plot = plot, width = 5.8, height = 10, units = "in", dpi = 200
)

