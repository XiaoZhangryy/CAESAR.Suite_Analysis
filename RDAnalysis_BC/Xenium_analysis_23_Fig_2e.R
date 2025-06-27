rm(list = ls())
################################################################################
# Figure 2.e
################################################################################

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

setwd(work_wd)

load("cols_cts_q.rda")

load_results <- function(sample, annotation_wd = annotation_wd) {
    result_file_name0 <- paste0(
        annotation_wd, gsub(" ", "_", sample), "_CAESAR.rda"
    )
    load(result_file_name0)
    anno.df_CAESAR <- anno.df
    
    result_file_name0 <- paste0(
        annotation_wd, gsub(" ", "_", sample), "_CelliD.rda"
    )
    load(result_file_name0)
    anno.df_CelliD <- anno.df

    result_file_name1 <- paste0(
        annotation_wd, gsub(" ", "_", sample),
        "_Seurat_scmap_SingleR_scPred.rda"
    )
    load(result_file_name1)
    anno.df_benchmark <- anno.df[, paste0(c(
        "Seurat", "Seuratunasg", "Seuratconf",
        "scmap", "scmapunasg", "scmapconf",
        "SingleR", "SingleRunasg", "SingleRconf",
        "scPred", "scPredunasg", "scPredconf"
    ), "_", gsub(" ", "_", sample))]

    anno.df <- cbind(anno.df_CAESAR, anno.df_benchmark)
    anno.df <- cbind(anno.df, anno.df_CelliD)
    return(anno.df)
}

sample <- "Sample 19"

methods <- c(
    "CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD"
)

plot_methods <- c(
  "iCAESAR", "iCAESARunasg", 
  paste0(methods, "_", gsub(" ", "_", sample)),
  paste0(methods, "unasg", "_", gsub(" ", "_", sample))
)

for (i in 1:4) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    annotation_wd <- paste0(ana_wd, "annotation_results/")
    setwd(ana_wd)
    
    load("seu.rda")

    coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
    if (i %in% c(1:2)) {
        plot_height <- coord[2] / 500
        plot_width <- coord[1] / 500
    } else {
        plot_height <- coord[2] / 1000
        plot_width <- coord[1] / 1000
    }

    anno.df <- load_results(sample, annotation_wd)
    seu <- AddMetaData(seu, anno.df, col.name = colnames(anno.df))

    for (pred in plot_methods) {
        Idents(seu) <- factor(seu[[pred]][, 1], levels = names(cols))
        plot <- DimPlot(seu, reduction = "pos", cols = cols, pt.size = 0.15, raster = FALSE) +
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

        plot2 <- DimPlot(
            seu,
            reduction = "pos", cols = cols, split.by = "ident",
            group.by = "ident", ncol = 3, pt.size = 0.2, raster = FALSE
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
            file = paste0(ana_wd, "DimPlot_pos_split_", pred, ".png"),
            plot = plot2, width = plot_width * 2, height = plot_height * 1.5, units = "in", dpi = 500
        )
    }


    plot3 <- FeaturePlot(
        seu,
        reduction = "pos", features = "iCAESARcf", pt.size = 0.15, raster = FALSE,
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
        file = paste0(ana_wd, "FeaturePlot_pos_", "iCAESARconf", ".png"),
        plot = plot3, width = plot_width, height = plot_height, units = "in", dpi = 500
    )

    message("section ", i, " done.")
}

