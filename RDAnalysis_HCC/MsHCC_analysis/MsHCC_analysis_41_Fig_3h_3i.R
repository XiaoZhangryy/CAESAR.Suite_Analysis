rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/MsHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

load("cols_cts_q.rda")
load("ms2hsgenes.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "seu.rda"))

    seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")
    y <- as.character(seu$CAESARua)
    stromal_cts <- c(
        "Fibroblasts", "Endothelial"
    )
    y[y %in% stromal_cts] <- "Stromal"
    immune_cts <- c(
        "B/Plasma", "Macrophages", "Neutrophil"
    )
    y[y %in% immune_cts] <- "Immu"
    y[y %in% c("Stromal", "Immu")] <- "Stromal&Immu"
    Idents(seu) <- y
    seu
})

cts <- c("HCC cell", "Stromal&Immu")

# ==============================================================================
# 1. calculate pvals
# ==============================================================================
load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

pathway_scoresList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "pathway_scores.rda"))
    pathway_scores
})

rgTestList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "df_rgTest.rda"))
    df_rgTest
})

enrich_pval <- apply(sapply(rgTestList, function(xx) {
    xx$wei.pval
}), 1, cauchy_combin, weight = sapply(seuList, ncol))
names(enrich_pval) <- pathway_names

enrich_pval_adj <- p.adjust(enrich_pval, method = "BH")

save(enrich_pval, file = paste0(enrich_wd, "enrich_pval.rda"))
save(enrich_pval_adj, file = paste0(enrich_wd, "enrich_pval_adj.rda"))

## celltype differential enriched
ct_enrich_pvalList <- lapply(1:4, function(i) {
    pathway_scores <- pathway_scoresList[[i]]
    y <- as.character(Idents(seuList[[i]]))
    ct_enrich(y, pathway_scores, cts)
})
save(ct_enrich_pvalList, file = paste0(enrich_wd, "ct_enrich_pvalList.rda"))

n_seu <- sapply(seuList, ncol)
ct_enrich_pval <- sapply(cts, function(ct) {
    n_seu_ct <- sapply(seuList, function(seu) {
        sum(Idents(seu) == ct)
    })
    pvals <- Reduce(cbind, lapply(ct_enrich_pvalList, function(xx) {
        xx[, ct]
    }))
    pvals <- apply(pvals, 1, cauchy_combin, weight = n_seu_ct)
    pvals
})
save(ct_enrich_pval, file = paste0(enrich_wd, "ct_enrich_pval.rda"))

## status differential enriched
pathway_scores <- Reduce(rbind, pathway_scoresList)
status <- c(
    rep("tumor", ncol(seuList[[1]])),
    rep("tumor", ncol(seuList[[2]])),
    rep("tumor-adj", ncol(seuList[[3]])),
    rep("tumor-adj", ncol(seuList[[4]]))
)
st_enrich_pval <- ct_enrich(status, pathway_scores, c("tumor", "tumor-adj"))
save(st_enrich_pval, file = paste0(enrich_wd, "st_enrich_pval.rda"))

# ==============================================================================
# 2. determine pathways
# ==============================================================================
load(paste0(enrich_wd, "enrich_pval_adj.rda"))
load(paste0(enrich_wd, "ct_enrich_pvalList.rda"))
load(paste0(enrich_wd, "ct_enrich_pval.rda"))
load(paste0(enrich_wd, "st_enrich_pval.rda"))

pathway_scoresList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "pathway_scores.rda"))
    pathway_scores
})

load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

gs <- names(enrich_pval_adj)[enrich_pval_adj < 0.05]
gs <- gs[!grepl("^GSE", gs)]
gs <- gs[sapply(gs, nchar) <= 80]

ct_enrich_pval_adj <- apply(ct_enrich_pval, 2, p.adjust, method = "BH")

# ------------------------------------------------------------------------------
gs <- names(enrich_pval_adj)[enrich_pval_adj < 0.05]

setNames(
    sapply(pathway_names_list, function(pp) length(intersect(pp, gs))),
    c("KEGG", "REACTOME", "CGP", "CM", "GOBP", "IMMUNESIGDB")
)

# ------------------------------------------------------------------------------
# pathway_used provided

ntop <- 6
pathway_used <- lapply(cts, function(ct, pathway_candidate = NULL) {
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
    pathway_candidate <- names(x)

    mean_scoreList <- lapply(seq_along(seuList), function(i) {
        pathway_scores <- pathway_scoresList[[i]]
        y <- as.character(Idents(seuList[[i]]))
        mean_score1 <- colMeans(pathway_scores[y == ct, pathway_candidate])
        mean_score2 <- colMeans(pathway_scores[y != ct, pathway_candidate])
        mean_score1 - mean_score2
    })
    mean_score <- Reduce(`+`, mapply(function(vec, weight) vec * weight, mean_scoreList, sapply(seuList, ncol), SIMPLIFY = FALSE))
    id_score <- order(mean_score, decreasing = TRUE)[seq_len(min(length(x), 100))]

    pathway_candidate <- pathway_candidate[id_score]
    x <- x[id_score]
    max_pval <- apply(sapply(ct_enrich_pvalList, function(pvalsi) {
        pvalsi[pathway_candidate, ct]
    }), 1, max)
    id <- order(max_pval, decreasing = FALSE)[seq_len(min(length(x), ntop))]

    return(x[id])
}, pathway_candidate = gs)
names(pathway_used) <- cts

pathway_usedList <- pathway_used

pathway_used <- pathway_usedList
cts_used <- cts
pathway_used <- lapply(pathway_used, function(p) names(p))
pathway_used <- unlist(pathway_used[cts_used], use.names = FALSE)
pathway_used <- unique(pathway_used)

save(pathway_used, file = paste0(enrich_wd, "pathway_used.rda"))

# ==============================================================================
# 3. plot pathways
# ==============================================================================
load(paste0(enrich_wd, "pathway_used.rda"))

## 3.1 data for plot
cts_pathway_scores <- Reduce(rbind, lapply(1:4, function(i) {
    y <- as.character(Idents(seuList[[i]]))
    pathway_scores <- pathway_scoresList[[i]][, pathway_used]
    cts_pathway_scores <- sapply(cts, function(ct) {
        colMeans(pathway_scores[y == ct, ])
    })
    rownames(cts_pathway_scores) <- pathway_used
    colnames(cts_pathway_scores) <- paste0(colnames(cts_pathway_scores), i)
    t(cts_pathway_scores)
}))
colnames(cts_pathway_scores) <- gsub("_", " ", colnames(cts_pathway_scores))

threshold <- 0.95
cts_pathway_pct <- Reduce(rbind, lapply(1:4, function(i) {
    y <- as.character(Idents(seuList[[i]]))
    pathway_scores <- pathway_scoresList[[i]][, pathway_used]
    cts_pathway_pct <- sapply(cts, function(ct) {
        colMeans(pathway_scores[y == ct, ] > threshold)
    })
    rownames(cts_pathway_pct) <- pathway_used
    colnames(cts_pathway_pct) <- paste0(colnames(cts_pathway_pct), i)
    t(cts_pathway_pct)
}))
colnames(cts_pathway_pct) <- gsub("_", " ", colnames(cts_pathway_pct))

reorder <- c(1, 3, 5, 7, 2, 4, 6, 8)
cts_pathway_scores <- cts_pathway_scores[reorder, ]
cts_pathway_pct <- cts_pathway_pct[reorder, ]
rownames(cts_pathway_scores) <- c(paste0("HCC cell(", 1:4, ")"), paste0("Stroma/Immune(", 1:4, ")"))
rownames(cts_pathway_pct) <- c(paste0("HCC cell(", 1:4, ")"), paste0("Stroma/Immune(", 1:4, ")"))

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
    j <- which(sapply(1:6, function(j) {
        pp0 %in% pathway_names_list[[j]]
    }))
    pp <- gsub("MODULE", "CANCER RELATED MODULE", pp)
    pp1 <- switch(j,
        paste0(pp),
        paste0(pp),
        paste0("CGP ", pp),
        paste0("CM ", pp),
        paste0(pp),
        paste0("IMMUNESIGDB", pp)
    )
    return(pp1)
})

levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)

if (FALSE) {
    plot <- data_melt %>%
        filter(Pct > 1) %>%
        ggplot(aes(x = Column, y = Row, color = Score, size = Pct)) +
        geom_point() +
        scale_size_continuous(
            breaks = c(10, 40, 70, 100), labels = c(10, 40, 70, 100),
            limits = c(0, 100)
        ) +
        # scale_color_viridis_c(name = "Average Score") +
        scale_color_gradientn(
            colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            values = scales::rescale(seq(0, 1, 0.25)),
            # colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            # values = scales::rescale(c(0, 0.3, 0.6, 0.9, 1.0)),
            limits = c(0, 1)
        ) +
        cowplot::theme_cowplot() +
        theme(
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12),
            # axis.title.x = element_text(margin = margin(t = 10), size = 14),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12, face = "bold"),
            # axis.title.y = element_text(margin = margin(r = 10), size = 14),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank(),
            legend.title = element_blank(),
            # legend.title = element_text(size = 12, face = "bold", margin = margin(b = 15)),
            # legend.title = element_blank(),
            legend.text = element_text(size = 14, face = "bold"),
            legend.spacing.y = unit(2.5, "cm"),
            legend.margin = margin(t = 2, unit = "cm")
        )

    ggsave(
        file = paste0(enrich_wd, "dotplot_cts_pathway_CAESARua_pval.png"),
        plot = plot, width = 8, height = 8, units = "in", dpi = 500
    )
} else {
    plot <- data_melt %>%
        filter(Pct > 1) %>%
        ggplot(aes(x = Row, y = Column, color = Score, size = Pct)) +
        geom_point() +
        scale_size_continuous(
            breaks = c(10, 40, 70, 100), labels = c(10, 40, 70, 100),
            limits = c(0, 100)
        ) +
        # scale_color_viridis_c(name = "Average Score") +
        scale_color_gradientn(
            colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            values = scales::rescale(seq(0, 1, 0.25)),
            # colors = c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#49006a"),
            # values = scales::rescale(c(0, 0.3, 0.6, 0.9, 1.0)),
            limits = c(0, 1)
        ) +
        cowplot::theme_cowplot() +
        theme(
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 12),
            # axis.title.x = element_text(margin = margin(t = 10), size = 14),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12, face = "bold"),
            # axis.title.y = element_text(margin = margin(r = 10), size = 14),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank(),
            legend.title = element_blank(),
            # legend.title = element_text(size = 12, face = "bold", margin = margin(b = 15)),
            # legend.title = element_blank(),
            legend.text = element_text(size = 14, face = "bold"),
            legend.spacing.y = unit(2.5, "cm"),
            legend.margin = margin(t = 1, unit = "cm")
        )

    ggsave(
        file = paste0(enrich_wd, "dotplot_cts_pathway_CAESARua_pval.png"),
        plot = plot, width = 8, height = 5, units = "in", dpi = 500
    )
}


## 3.3 spatial plot of each pathway on each section
features <- pathway_used
for (i in 1:4) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    seu <- seuList[[i]]

    coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
    plot_height <- coord[2] / 500
    plot_width <- coord[1] / 500

    seu <- AddMetaData(seu, as.data.frame(pathway_scoresList[[i]]))

    for (i_feature in seq_along(features)) {
        feature <- features[i_feature]

        plot <- FeaturePlot(seu, features = feature, reduction = "pos", pt.size = 1.9) +
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
            file = paste0(ana_wd, "FeaturePlot_pos_", feature, "_HCC", i, "_pure.png"),
            plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
        )
    }


    for (jj in seq_along(cts)) {
        ct <- cts[jj]
        kegg_wd_i_ct <- paste0(ana_wd, gsub("/|&| ", "_", ct), "/")
        dir.create(kegg_wd_i_ct, showWarnings = FALSE)

        features_ct <- pathway_used[1:6 + (jj - 1) * 6]
        for (feature in features_ct) {
            file <- paste0(ana_wd, "FeaturePlot_pos_", feature, "_HCC", i, "_pure.png")
            file0 <- paste0(kegg_wd_i_ct, "FeaturePlot_pos_", feature, "_HCC", i, "_pure.png")
            system(paste0("mv ", file, " ", file0))
        }
    }
}



## 3.5 top pathways of each section
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtree)

rgTestList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "df_rgTest.rda"))
    df_rgTest
})

for (i in 1:4) {
    ntop <- 5

    enrich_pval_adj <- p.adjust(rgTestList[[i]]$wei.pval, method = "BH")
    gs <- pathway_names[enrich_pval_adj < 0.05]
    gs <- gs[!grepl("^GSE", gs)]
    gs <- gs[sapply(gs, nchar) <= 80]

    ntop <- 6
    pathway_used <- lapply(cts, function(ct, pathway_candidate = NULL) {
        if (is.null(pathway_candidate)) {
            pathway_candidate <- rownames(ct_enrich_pval_adj)
        }

        ct_enrich_pval_adj <- p.adjust(ct_enrich_pvalList[[i]][, ct], method = "BH")

        x <- ct_enrich_pval_adj[pathway_candidate]
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
        # when p values are the same, sort by difference and max of p values
        pathway_candidate <- names(x)

        pathway_scores <- pathway_scoresList[[i]]
        y <- as.character(Idents(seuList[[i]]))
        mean_score1 <- colMeans(pathway_scores[y == ct, pathway_candidate])
        mean_score2 <- colMeans(pathway_scores[y != ct, pathway_candidate])
        mean_score <- mean_score1 - mean_score2
        id_score <- order(mean_score, decreasing = TRUE)[seq_len(min(length(x), ntop))]

        return(x[id_score])
    }, pathway_candidate = gs)
    names(pathway_used) <- cts

    pathway_usedList <- pathway_used

    cts_used <- cts
    pathway_used <- lapply(pathway_used, function(p) names(p))
    pathway_used <- unlist(pathway_used[cts_used], use.names = FALSE)
    pathway_used <- unique(pathway_used)

    ### cts_pathway_scores
    y <- as.character(Idents(seuList[[i]]))
    pathway_scores <- pathway_scoresList[[i]][, pathway_used]

    cts_pathway_scores <- t(sapply(cts, function(ct) {
        colMeans(pathway_scores[y == ct, ])
    }))
    colnames(cts_pathway_scores) <- gsub("_", " ", pathway_used)
    rownames(cts_pathway_scores) <- c("HCC cell", "Stroma/Immune")

    ### cts_pathway_pct
    threshold <- 0.95
    cts_pathway_pct <- t(sapply(cts, function(ct) {
        colMeans(pathway_scores[y == ct, ] > threshold)
    }))
    colnames(cts_pathway_pct) <- gsub("_", " ", pathway_used)
    rownames(cts_pathway_pct) <- c("HCC cell", "Stroma/Immune")

    ### data_melt
    data_score <- reshape2::melt(cts_pathway_scores)
    colnames(data_score) <- c("Row", "Column", "Score")
    data_pct <- reshape2::melt(cts_pathway_pct * 100)
    colnames(data_pct) <- c("Row", "Column", "Pct")

    data_melt <- left_join(data_score, data_pct, by = c("Row", "Column"))
    data_melt$Row <- factor(data_melt$Row, levels = rownames(cts_pathway_scores))
    data_melt$Column <- factor(data_melt$Column, levels = colnames(cts_pathway_scores))


    save(pathway_usedList, cts_pathway_scores, cts_pathway_pct, data_melt, file = paste0(enrich_wd, "plot_data_List_pval_section", i, ".rda"))

    ### dotplot
    levels(data_melt$Column) <- sapply(colnames(cts_pathway_scores), function(pp) {
        pp0 <- gsub(" ", "_", pp)
        j <- which(sapply(1:6, function(j) {
            pp0 %in% pathway_names_list[[j]]
        }))
        pp <- gsub("MODULE", "CANCER RELATED MODULE", pp)
        pp1 <- switch(j,
            paste0(pp),
            paste0(pp),
            paste0("CGP ", pp),
            paste0("CM ", pp),
            paste0(pp),
            paste0("IMMUNESIGDB", pp)
        )
        return(pp1)
    })

    levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)

    colnames(data_melt) <- c("Row", "Column", "Ave. Score", "Pct. Enriched")

    plot <- data_melt %>%
        filter(`Pct. Enriched` > 1) %>%
        ggplot(aes(x = Column, y = Row, color = `Ave. Score`, size = `Pct. Enriched`)) +
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
            legend.title = element_text(size = 16, face = "bold"),
            # legend.title = element_text(size = 12, face = "bold", margin = margin(b = 15)),
            # legend.title = element_blank(),
            legend.text = element_text(size = 14, face = "bold"),
            legend.spacing.y = unit(-5, "cm"),
            legend.margin = margin(t = 6, unit = "cm")
        )
    # labs(x = "KEGG pathways", y = "SAX predictions")

    ggsave(
        file = paste0(enrich_wd, "dotplot_cts_section_pathway_matrix_CAESARua_pval_HCC", i, ".png"),
        plot = plot, width = 7, height = 5, units = "in", dpi = 200
    )
}
