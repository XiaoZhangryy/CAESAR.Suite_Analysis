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
# 1. spatial plot - fig 2.e
################################################################################
for (i in 1:4) {
    seu <- seuList[[i]]
    ana_wd <- paste0(work_wd, "BC", i, "/")

    coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
    if (i %in% c(1:2)) {
        plot_height <- coord[2] / 500
        plot_width <- coord[1] / 500
    } else {
        plot_height <- coord[2] / 1000
        plot_width <- coord[1] / 1000
    }

    for (pred in c("iCAESARunasg", "CelliDunasg_Sample_19")) {
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
        file = paste0(ana_wd, "FeaturePlot_pos_", "iCAESARcf", ".png"),
        plot = plot3, width = plot_width, height = plot_height, units = "in", dpi = 500
    )

    message("section ", i, " done.")
}

################################################################################
# 1. co-embedding plot - fig. 2.c
################################################################################
load("markerList.rda")
load("sgInt.rda")
load("seu_scList.rda")
load("sg_List.rda")
load("ct_ratio.rda")
load("sg_sc_List.rda")

marker_freq <- setNames(lapply(celltypes, function(ct) {
    genes <- lapply(markerList, function(marker) marker[[ct]]) %>%
        unlist() %>%
        table()
    # sort(decreasing = TRUE)
    id <- sapply(names(genes), function(gene) match(gene, sgInt[[ct]]))
    names(genes)[order(-genes, id, decreasing = FALSE)]
}), celltypes)

for (i in 1:4) {
    seu <- seuList[[i]]
    Idents(seu) <- seu$iCAESARunasg
    ana_wd <- paste0(work_wd, "BC", i, "/")

    ct_used <- names(ct_ratio[[i]][ct_ratio[[i]] > 5e-3])
    ct_used <- intersect(ct_used, celltypes)

    df_coumap_gene <- Reduce(rbind, lapply(
        ct_used, function(ct) {
            sig_sp <- sg_List[[i]][[ct]]
            gene_used <- marker_freq[[ct]][1:2]
            sig_sp[match(gene_used, sig_sp$gene), ]
        }
    ))
    coumap_gene <- df_coumap_gene %>%
        .$gene %>%
        unique()

    seu <- CoUMAP(
        seu,
        reduction = "caesar", reduction.name = "CoUMAP",
        gene.set = coumap_gene, n_neighbors = 50
    )

    plot3 <- CoUMAP.plot(
        seu,
        reduction = "CoUMAP", gene_txtdata = df_coumap_gene,
        cols = c("gene" = "#000000", cols),
        pt_size = 0.2, pt_text_size = 4, alpha = 0.8
    ) + guides(
        shape = guide_legend(override.aes = list(size = 3)),
        color = guide_legend(ncol = 1, override.aes = list(size = 3))
    )
    ggsave(
        file = paste0(ana_wd, "CoUMAP_iCAESARunasg_BC", i, "_marker.png"),
        plot = plot3, width = 9, height = 5, units = "in", dpi = 500
    )
}

figs_wd_marker <- paste0(figs_wd, "CoUMAP_marker/")
dir.create(figs_wd_marker, showWarnings = FALSE)
for (i in seq_along(seu_scList)) {
    seu_sc <- seu_scList[[i]]
    sample <- gsub(" ", "_", names(seu_scList)[[i]])

    sc_ct_ratio <- table(seu_sc$merged_clusters)
    sc_ct_ratio <- sc_ct_ratio / sum(sc_ct_ratio)

    ct_used <- names(sc_ct_ratio[sc_ct_ratio > 5e-3])
    ct_used <- intersect(ct_used, celltypes)

    df_coumap_gene_sc <- Reduce(rbind, lapply(
        ct_used, function(ct) {
            gene <- intersect(marker_freq[[ct]], markerList[[i]][[ct]])[1:2]
            sig_sp <- sg_sc_List[[i]][[ct]]
            sig_sp[match(gene, sig_sp$gene), ]
        }
    ))
    coumap_gene_sc <- df_coumap_gene_sc %>%
        .$gene %>%
        unique()

    Idents(seu_sc) <- factor(seu_sc$merged_clusters, levels = celltypes)

    set.seed(1)
    seu_sc <- CoUMAP(
        seu_sc,
        reduction = "ncfm", reduction.name = "CoUMAP",
        gene.set = coumap_gene_sc, n_neighbors = 30
    )

    plot3 <- CoUMAP.plot(
        seu_sc,
        reduction = "CoUMAP", gene_txtdata = df_coumap_gene_sc,
        cols = c("gene" = "#000000", cols),
        pt_size = 0.6, pt_text_size = 4, alpha = 0.8
    ) + guides(
        shape = guide_legend(override.aes = list(size = 3)),
        color = guide_legend(ncol = 1, override.aes = list(size = 3))
    )
    ggsave(
        file = paste0(figs_wd_marker, "sc_CoUMAP_marker_freq_", sample, ".png"),
        plot = plot3, width = 9, height = 5, units = "in", dpi = 500
    )

    plot5 <- DimPlot(
        seu_sc,
        reduction = "CoUMAP", cols = cols, group.by = "ident", raster = FALSE
    ) +
        NoLegend() +
        theme(
            axis.line = element_blank(), # Remove the axis lines
            axis.text.x = element_blank(), # Remove the text on the x-axis
            axis.text.y = element_blank(), # Remove the text on the y-axis
            axis.ticks = element_blank(), # Remove the ticks
            axis.title.x = element_blank(), # Remove the x-axis title
            axis.title.y = element_blank(),
            plot.title = element_blank()
        ) +
        annotate("text",
            x = Inf, y = -Inf, label = names(seu_scList)[[i]],
            hjust = 1.1, vjust = -0.5, size = 14, color = "black"
        )
    ggsave(
        file = paste0(figs_wd_marker, "sc_UMAP_marker_freq_", sample, ".png"),
        plot = plot5, width = 9, height = 5, units = "in", dpi = 200
    )
}

