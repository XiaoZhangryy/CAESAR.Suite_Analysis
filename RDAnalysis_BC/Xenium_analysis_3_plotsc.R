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
load("seu_scList.rda")
load("seu_sc.rda")

################################################################################
# 1. fig 2a - umap plot of scRNA-seq
################################################################################
seu_sc <- RunUMAP(
    seu_sc,
    reduction = "ncfm", dims = 1:q_est, reduction.name = "umap_ncfm"
)

match_sample <- setNames(names(seu_scList), sapply(seu_scList, function(seu) seu$orig.ident[1]))
seu_sc$samplenumber <- factor(match_sample[seu_sc$orig.ident], levels = match_sample)

plot1 <- DimPlot(
    seu_sc,
    reduction = "umap_ncfm", group.by = "samplenumber",
    raster = FALSE, pt.size = 0.01,
    cells = sample(colnames(seu_sc), 20000),
    shuffle = TRUE
) +
    theme(
        plot.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(face = "bold", size = 16),
        plot.background = element_blank()
    ) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))
ggsave(
    file = paste0(figs_wd, "sc_DimPlot_umap_ncfm.png"),
    plot = plot1, width = 10, height = 4, units = "in", dpi = 500
)

################################################################################
# 2. fig 2b - Stacked column chart
################################################################################
sc_ct_ratios <- sapply(seu_scList, function(seu) {
    counts <- table(seu$merged_clusters)
    ratio <- setNames(rep(0, length(celltypes)), celltypes)
    ratio[names(counts)] <- counts / sum(counts)
    ratio
})
sc_ct_ratios <- sc_ct_ratios[c(
    "Cancer Epithelial", "T-cells", "Myeloid", "Endothelial", "CAFs", "PVL",
    "B-cells", "Plasmablasts", "Normal Epithelial"
), ]
cols_scc <- cols[rownames(sc_ct_ratios)]
# rownames(sc_ct_ratios) <- c(
#     "CE", "T-cells", "Myeloid", "Endothelial", "CAFs", "PVL",
#     "B-cells", "Plasmablasts", "NE"
# )
# names(cols_scc) <- rownames(sc_ct_ratios)

# Convert matrix to data frame
data_df <- as.data.frame(t(sc_ct_ratios))
data_df$Sample <- rownames(data_df)

df_cn <- colnames(data_df)
df_cn[df_cn == "Cancer Epithelial"] <- "Cancer epithelial"
df_cn[df_cn == "Normal Epithelial"] <- "Normal epithelial"
colnames(data_df) <- df_cn

# Reshape data frame
data_melted <- reshape2::melt(data_df, id.vars = "Sample", variable.name = "Cell type", value.name = "Ratio")
data_melted$Sample <- factor(data_melted$Sample, levels = rownames(data_df))
# data_melted$`Cell type` <- factor(data_melted$`Cell type`, levels = rev(celltypes[c(1, 9, 7, 5, 2, 4, 8, 6, 3)]))
data_melted$`Cell type` <- factor(data_melted$`Cell type`, levels = c(
    "Cancer epithelial", "Normal epithelial", "PVL",
    "Endothelial", "CAFs", "Myeloid", "Plasmablasts", "T-cells", "B-cells"
))

cols_scc <- cols_scc[rev(celltypes[c(1, 9, 7, 5, 2, 4, 8, 6, 3)])]
names(cols_scc) <- levels(data_melted$`Cell type`)

# Create the stacked column chart
plot2 <- ggplot(data_melted, aes(x = Sample, y = Ratio, fill = `Cell type`)) +
    geom_bar(stat = "identity") + # Create stacked bar chart
    scale_fill_manual(values = cols_scc) +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # Background with border
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = "bottom", # Position legend at the bottom
        legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15)
    ) +
    labs(title = NULL, x = NULL, y = NULL)

ggsave(
    file = paste0(figs_wd, "sc_Stacked_Column_Chart.png"),
    plot = plot2, width = 10, height = 5, units = "in", dpi = 200
)














