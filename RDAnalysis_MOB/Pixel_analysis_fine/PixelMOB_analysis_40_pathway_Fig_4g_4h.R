rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")

load(paste0(ana_wd, "seu.rda"))
load(paste0(ana_wd, "pathway_scores.rda"))
load(paste0(ana_wd, "df_rgTest.rda"))

y <- as.character(seu$CAESARunasg)
cts <- c("GC", "M/TC", "PGC", "Neuron.OSN", "OEC", "Mes")

# ==============================================================================
# 1. calculate pvals
# ==============================================================================
load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

ct_enrich_pval <- ct_enrich(y, pathway_scores, cts)

save(ct_enrich_pval, file = paste0(enrich_wd, "ct_enrich_pval.rda"))

enrich_pval <- df_rgTest$wei.pval
enrich_pval_adj <- p.adjust(enrich_pval, method = "BH")
names(enrich_pval_adj) <- names(enrich_pval) <- pathway_names

save(enrich_pval, file = paste0(enrich_wd, "enrich_pval.rda"))
save(enrich_pval_adj, file = paste0(enrich_wd, "enrich_pval_adj.rda"))

# ==============================================================================
# 2. determine pathways
# ==============================================================================
# load(paste0(enrich_wd, "enrich_pval.rda"))
load(paste0(enrich_wd, "enrich_pval_adj.rda"))
load(paste0(enrich_wd, "ct_enrich_pval.rda"))

load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

# gs <- colnames(pathway_scoresList[[1]])
gs <- names(enrich_pval_adj)[enrich_pval_adj < 0.05]
gs <- gs[!grepl("^GSE", gs)]
gs <- gs[sapply(gs, nchar) <= 80]

ct_enrich_pval_adj <- apply(ct_enrich_pval, 2, p.adjust, method = "BH")


ntop <- 5
pathway_used <- lapply(cts, function(ct, y, pathway_scores, pathway_candidate = NULL) {
    if (is.null(pathway_candidate)) {
        pathway_candidate <- rownames(ct_enrich_pval_adj)
    }

    x <- ct_enrich_pval_adj[pathway_candidate, ct]
    x <- sort(x, decreasing = FALSE)
    x <- x[x <= 0.05]
    if (length(x) == 0) {
        return(NULL)
    }
    pathway <- x[seq_len(min(length(x), ntop))]
    if (max(pathway) != 0) {
        return(pathway)
    }
    x <- x[x == 0]
    # when p values are the same, sort by difference
    pathway_candidate <- names(x)

    mean_score <- colMeans(pathway_scores[y == ct, pathway_candidate]) - colMeans(pathway_scores[y != ct, pathway_candidate])
    id_score <- order(mean_score, decreasing = TRUE)
    # id_score <- id_score[seq_len(min(length(x), ntop))]

    return(names(x[id_score]))
}, pathway_candidate = gs, y = y, pathway_scores = pathway_scores)
names(pathway_used) <- cts

pathway_used[["GC"]] <- pathway_used[["GC"]][c(5, 11, 14, 1, 58)]
pathway_used[["M/TC"]] <- pathway_used[["M/TC"]][c(3, 5, 6, 11, 25)]
pathway_used[["PGC"]] <- pathway_used[["PGC"]][c(6, 7, 8, 21, 30)]
pathway_used[["Neuron.OSN"]] <- pathway_used[["Neuron.OSN"]][c(68, 9, 1, 4, 92)]
pathway_used[["OEC"]] <- pathway_used[["OEC"]][c(2, 11, 29, 33, 20)]
pathway_used[["Mes"]] <- pathway_used[["Mes"]][c(21, 24, 82, 117, 32)]



pathway_usedList <- pathway_used

cts_used <- rev(cts)
# pathway_used <- lapply(pathway_used, function(p) names(p))
pathway_used <- unlist(pathway_usedList[cts_used], use.names = FALSE)
pathway_used <- unique(pathway_used)

save(pathway_used, file = paste0(enrich_wd, "pathway_used.rda"))


# ==============================================================================
# 3. plot pathways
# ==============================================================================
load(paste0(enrich_wd, "pathway_used.rda"))
## 3.1 data for plot
cts_pathway_scores <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_used])
}))
colnames(cts_pathway_scores) <- gsub("_", " ", colnames(cts_pathway_scores))
rownames(cts_pathway_scores) <- c("Mes", "OEC", "OSNs", "PGC", "M/TC", "GC")

threshold <- 0.95
cts_pathway_pct <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_used] > threshold)
}))
colnames(cts_pathway_pct) <- gsub("_", " ", colnames(cts_pathway_pct))
rownames(cts_pathway_pct) <- c("Mes", "OEC", "OSNs", "PGC", "M/TC", "GC")

data_score <- reshape2::melt(cts_pathway_scores)
colnames(data_score) <- c("Row", "Column", "Score")
data_pct <- reshape2::melt(cts_pathway_pct * 100)
colnames(data_pct) <- c("Row", "Column", "Pct")

data_melt <- left_join(data_score, data_pct, by = c("Row", "Column"))
data_melt$Row <- factor(data_melt$Row, levels = rownames(cts_pathway_scores))
data_melt$Column <- factor(data_melt$Column, levels = colnames(cts_pathway_scores))

save(cts_pathway_scores, cts_pathway_pct, data_melt, file = paste0(enrich_wd, "plot_data_pval.rda"))


## 3.2 pathway dotplot
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtree)

levels(data_melt$Column) <- sapply(colnames(cts_pathway_scores), function(pp) {
    pp0 <- gsub(" ", "_", pp)
    j <- which(sapply(seq_along(pathway_names_list), function(j) {
        pp0 %in% pathway_names_list[[j]]
    }))
    pp <- gsub("MODULE", "CANCER RELATED MODULE", pp)
    pp1 <- switch(j,
        paste0("H ", pp),
        paste0(pp),
        paste0(pp),
        paste0("GTRD ", pp),
        paste0(pp),
        paste0(pp),
        paste0(pp)
    )
    return(pp1)
})

levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)

plot <- data_melt %>%
    filter(Pct > 1) %>%
    ggplot(aes(x = Column, y = Row, color = Score, size = Pct)) +
    geom_point() +
    # scale_color_viridis_c(name = "Average Score") +
    scale_size_continuous(
        breaks = c(10, 40, 70, 100), labels = c(10, 40, 70, 100),
        limits = c(0, 100)
    ) +
    scale_color_gradientn(
        colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
        values = scales::rescale(seq(0, 1, 0.25)),
        limits = c(0, 1)
    ) +
    cowplot::theme_cowplot() +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
        # axis.title.x = element_text(margin = margin(t = 10), size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        # axis.title.y = element_text(margin = margin(r = 10), size = 14),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        # legend.title = element_text(size = 12, face = "bold", margin = margin(b = 15)),
        # legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "bold"),
        legend.spacing.y = unit(1.5, "cm"),
        legend.margin = margin(t = 2.7, unit = "cm")
    )
# labs(x = "KEGG pathways", y = "SAX predictions")

ggsave(
    file = paste0(enrich_wd, "dotplot_cts_pathway_matrix_CAESARunasg_pval.png"),
    plot = plot, width = 15, height = 7, units = "in", dpi = 200
)



## 3.3 spatial plot of each pathway on each section
features <- pathway_used

seu <- AddMetaData(seu, as.data.frame(pathway_scores))

coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 2000
plot_width <- coord[1] / 2000

for (i_feature in seq_along(features)) {
    feature <- features[i_feature]

    plot <- FeaturePlot(seu, features = feature, reduction = "pos", raster = FALSE) +
        scale_color_gradientn(
            # colors = c("#f6eff7", "#feebe2", "#f768a1", "#7a0177", "#6e016b"),
            colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            values = scales::rescale(seq(0, 1, 0.25)),
            limits = c(0, 1)
        ) +
        theme(
            legend.position = "right",
            legend.justification = "center",
            legend.box = "vertical"
        )
    ggsave(
        file = paste0(enrich_wd, "FeaturePlot_pos_", feature, ".png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
    )

    plot <- FeaturePlot(seu, features = feature, reduction = "pos", raster = FALSE) +
        scale_color_gradientn(
            # colors = c("#f6eff7", "#feebe2", "#f768a1", "#7a0177", "#6e016b"),
            colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            values = scales::rescale(seq(0, 1, 0.25)),
            limits = c(0, 1)
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none"
            # legend.justification = "center",
            # legend.box = "vertical"
        ) + # Position legend at the bottom
        labs(title = levels(data_melt$Column)[i_feature], x = NULL, y = NULL)

    ggsave(
        file = paste0(enrich_wd, "FeaturePlot_pos_", feature, "_pure.png"),
        plot = plot, width = plot_width, height = plot_height * 1.1, units = "in", dpi = 200
    )
}


for (i_feature in seq_along(features)) {
    feature <- features[i_feature]

    plot <- FeaturePlot(seu, features = feature, reduction = "pos", raster = FALSE) +
        scale_color_gradientn(
            # colors = c("#f6eff7", "#feebe2", "#f768a1", "#7a0177", "#6e016b"),
            colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            values = scales::rescale(seq(0, 1, 0.25)),
            limits = c(0, 1)
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none"
            # legend.justification = "center",
            # legend.box = "vertical"
        ) + # Position legend at the bottom
        labs(title = NULL, x = NULL, y = NULL)

    ggsave(
        file = paste0(enrich_wd, "FeaturePlot_pos_", feature, "_none_title.png"),
        plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
    )
}



for (ct in cts) {
    enrich_wd_ct <- paste0(enrich_wd, gsub("-|/", "_", ct), "/")
    dir.create(enrich_wd_ct, showWarnings = FALSE)

    features_ct <- pathway_usedList[[ct]]
    for (feature in features_ct) {
        file <- paste0(enrich_wd, "FeaturePlot_pos_", feature, "_pure.png")
        file0 <- paste0(enrich_wd_ct, "FeaturePlot_pos_", feature, "_pure.png")
        system(paste0("cp ", file, " ", file0))
    }
}



