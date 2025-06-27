rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")
load(paste0(ana_wd, "seu.rda"))

################################################################################
# 1. co-embedding plot w.r.t. markerlist
################################################################################
load(paste0(ana_wd, "sg_list.rda"))
load("markerList.rda")
load("seu_sc.rda")
load("sg_sc_List.rda")

cts_ratio <- table(seu$CAESARunasg) / ncol(seu)
cts_used <- names(cts_ratio)[cts_ratio > 1e-2]

marker <- setNames(lapply(cts_used, function(ct) {
    genes <- markerList[[1]][[ct]]
    id <- sapply(genes, function(gene) {
        match(gene, sg_list[[ct]]$gene)
    })
    genes[order(id, decreasing = FALSE)][1:2]
}), cts_used)



Idents(seu) <- seu$CAESARunasg
df_coumap_gene <- Reduce(rbind, lapply(
    cts_used, function(ct) {
        sig_sp <- sg_list[[ct]]
        sig_sp[match(marker[[ct]], sig_sp$gene), ]
    }
))
coumap_gene <- df_coumap_gene %>%
    .$gene %>%
    unique()

seu <- CoUMAP(
    seu,
    reduction = "caesar", reduction.name = "CoUMAP",
    gene.set = coumap_gene, n_neighbors = 30
)

plot3 <- CoUMAP.plot(
    seu,
    reduction = "CoUMAP", gene_txtdata = df_coumap_gene,
    cols = c("gene" = "#000000", cols),
    pt_size = 0.4, pt_text_size = 2, alpha = 0.8
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
)
ggsave(
    file = paste0(ana_wd, "CoUMAP_CAESARunasg_marker.png"),
    plot = plot3, width = 9, height = 5, units = "in", dpi = 200
)

plot4 <- CoUMAP.plot(
    seu,
    reduction = "CoUMAP",
    cols = c("gene" = "#000000", cols),
    pt_size = 0.4, pt_text_size = 2, alpha = 0.8
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
) + mytheme()
ggsave(
    file = paste0(ana_wd, "CoUMAP_CAESARunasg_marker_blank.png"),
    plot = plot4, width = 4.2, height = 3.5, units = "in", dpi = 200
)




df_coumap_gene_sc <- Reduce(rbind, lapply(
    cts_used, function(ct) {
        sig_sc <- sg_sc_List[[1]][[ct]]
        sig_sc[match(marker[[ct]], sig_sc$gene), ]
    }
))

coumap_gene_sc <- df_coumap_gene_sc %>%
    .$gene %>%
    unique()

Idents(seu_sc) <- factor(seu_sc$merged_clusters, levels = names(cols))

seu_sc <- CoUMAP(
    seu_sc,
    reduction = "ncfm", reduction.name = "CoUMAP",
    gene.set = coumap_gene_sc, n_neighbors = 30, min_dist = 0.3, spread = 1
)

plot3 <- CoUMAP.plot(
    seu_sc,
    reduction = "CoUMAP", gene_txtdata = df_coumap_gene_sc,
    cols = c("gene" = "#000000", cols),
    pt_size = 0.2, pt_text_size = 2, alpha = 0.8
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
)
ggsave(
    file = paste0(ana_wd, "sc_CoUMAP_CAESARunasg_marker.png"),
    plot = plot3, width = 9, height = 5, units = "in", dpi = 200
)

plot4 <- CoUMAP.plot(
    seu_sc,
    reduction = "CoUMAP",
    cols = c("gene" = "#000000", cols),
    pt_size = 0.2, pt_text_size = 2, alpha = 0.8
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
) + mytheme()
ggsave(
    file = paste0(ana_wd, "sc_CoUMAP_CAESARunasg_marker_blank.png"),
    plot = plot4, width = 4.2, height = 3.5, units = "in", dpi = 200
)


################################################################################
# 2. plot marker genes
################################################################################
coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 2000
plot_width <- coord[1] / 2000

genes <- coumap_gene
for (gene in genes) {
    plot <- FeaturePlot(seu, features = gene, reduction = "pos", min.cutoff = "q10", max.cutoff = "q80", cols = c("#fffedf", "#d7301f"), pt.size = 0.4, order = TRUE) +
        NoLegend() +
        theme(
            axis.line = element_blank(), # Remove the axis lines
            axis.text.x = element_blank(), # Remove the text on the x-axis
            axis.text.y = element_blank(), # Remove the text on the y-axis
            axis.ticks = element_blank(), # Remove the ticks
            axis.title.x = element_blank(), # Remove the x-axis title
            axis.title.y = element_blank(),
            plot.title = element_blank()
        )
    ggsave(
        file = paste0(figs_wd, "FeaturePlot_pos_", gene, ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
    )
}

