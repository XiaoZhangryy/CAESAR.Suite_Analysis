rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)

source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

load("cols_cts_q.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    Idents(seu) <- factor(seu$iCAESARunasg, levels = celltypes)
    seu
})


ctsList <- lapply(1:4, function(i) {
    seu <- seuList[[i]]
    cts.ratio <- table(Idents(seu)) / ncol(seu)
    cts <- setdiff(names(cts.ratio[cts.ratio > 5e-3]), "unassigned")
    cts
})
cts <- celltypes[c(1, 9, 7, 5, 2, 4, 8, 6, 3)]

# ==============================================================================
# 1. calculate pvals
# ==============================================================================
load("pathwayList.rda")

pathway_scoresList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "pathway_scores.rda"))
    pathway_scores
})

pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

ct_enrich_pvalList <- lapply(1:4, function(i) {
    pathway_scores <- pathway_scoresList[[i]]
    CAESAR.CTDEP(seuList[[i]], ident = "iCAESARunasg", pathway_scores, cts = cts)
})
save(ct_enrich_pvalList, file = paste0(enrich_wd, "ct_enrich_pvalList.rda"))


n_seu <- sapply(seuList, ncol)
ct_enrich_pval <- sapply(cts, function(ct) {
    n_seu_ct <- sapply(seuList, function(seu) {
        sum(seu$iCAESARunasg == ct)
    })
    pvals <- Reduce(cbind, lapply(ct_enrich_pvalList, function(xx) {
        xx[, ct]
    }))
    pvals <- apply(pvals, 1, Cauchy.Combination, weight = n_seu_ct)
    pvals
})
save(ct_enrich_pval, file = paste0(enrich_wd, "ct_enrich_pval.rda"))


rgTestList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "df_rgTest.rda"))
    df_rgTest
})

enrich_pval <- apply(sapply(rgTestList, function(xx) {
    xx$wei.pval
}), 1, Cauchy.Combination, weight = sapply(seuList, ncol))
names(enrich_pval) <- pathway_names

enrich_pval_adj <- p.adjust(enrich_pval, method = "BH")

save(enrich_pval, file = paste0(enrich_wd, "enrich_pval.rda"))
save(enrich_pval_adj, file = paste0(enrich_wd, "enrich_pval_adj.rda"))

# ==============================================================================
# 2. determine pathways
# ==============================================================================
load(paste0(enrich_wd, "enrich_pval_adj.rda"))
load(paste0(enrich_wd, "ct_enrich_pvalList.rda"))
load(paste0(enrich_wd, "ct_enrich_pval.rda"))

pathway_scoresList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "pathway_scores.rda"))
    pathway_scores
})

load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

# gs <- colnames(pathway_scoresList[[1]])
gs <- names(enrich_pval_adj)[enrich_pval_adj < 0.05]
gs <- gs[!grepl("^GSE", gs)]
gs <- gs[sapply(gs, nchar) <= 80]

ct_enrich_pval_adj <- apply(ct_enrich_pval, 2, p.adjust, method = "BH")


ntop <- 5
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
    return(pathway)
}, pathway_candidate = gs)
names(pathway_used) <- cts


pathway_usedList <- pathway_used
cts_used <- cts
pathway_used <- lapply(pathway_used, function(p) names(p))
pathway_used <- unlist(pathway_used[cts_used], use.names = FALSE)
pathway_used <- unique(pathway_used)

save(pathway_used, file = paste0(enrich_wd, "pathway_used.rda"))


# ==============================================================================
# 3. plot pathways - fig. 2.g and 2.h
# ==============================================================================
load(paste0(enrich_wd, "pathway_used.rda"))

## 3.1 data for plot
n_cts <- sapply(cts, function(ct) {
    sum(sapply(1:4, function(i) {
        if (ct %in% ctsList[[i]]) {
            return(sum(seuList[[i]]$iCAESARunasg == ct))
        } else {
            return(0)
        }
    }))
})

cts_pathway_scoresList <- lapply(1:4, function(i) {
    cts0 <- ctsList[[i]][ctsList[[i]] %in% cts]
    y <- seuList[[i]]$iCAESARunasg
    pathway_scores <- pathway_scoresList[[i]]
    cts_pathway_scores <- matrix(0, length(cts), length(pathway_used), dimnames = list(
        cts, pathway_used
    ))
    cts_pathway_scores[cts0, ] <- t(sapply(cts0, function(ct) {
        colSums(pathway_scores[y == ct, pathway_used])
    }))
    cts_pathway_scores
})
cts_pathway_scores <- Reduce(`+`, cts_pathway_scoresList) / n_cts
colnames(cts_pathway_scores) <- gsub("_", " ", colnames(cts_pathway_scores))

threshold <- 0.95
cts_pathway_pctList <- lapply(1:4, function(i) {
    cts0 <- ctsList[[i]][ctsList[[i]] %in% cts]
    y <- seuList[[i]]$iCAESARunasg
    pathway_scores <- pathway_scoresList[[i]]
    cts_pathway_pct <- matrix(0, length(cts), length(pathway_used), dimnames = list(
        cts, pathway_used
    ))
    cts_pathway_pct[cts0, ] <- t(sapply(cts0, function(ct) {
        colSums(pathway_scores[y == ct, pathway_used] > threshold)
    }))
    cts_pathway_pct
})
cts_pathway_pct <- Reduce(`+`, cts_pathway_pctList) / n_cts
colnames(cts_pathway_pct) <- gsub("_", " ", colnames(cts_pathway_pct))


data_score <- reshape2::melt(cts_pathway_scores)
colnames(data_score) <- c("Row", "Column", "Score")
data_pct <- reshape2::melt(cts_pathway_pct * 100)
colnames(data_pct) <- c("Row", "Column", "Pct")

data_melt <- left_join(data_score, data_pct, by = c("Row", "Column"))
data_melt$Row <- factor(data_melt$Row, levels = rownames(cts_pathway_scores))
data_melt$Column <- factor(data_melt$Column, levels = colnames(cts_pathway_scores))

save(cts_pathway_scores, cts_pathway_pct, data_melt, file = paste0(enrich_wd, "plot_data_pval.rda"))

replace_cts <- function(cts) {
    lv_cts <- levels(cts)
    lv_cts[lv_cts == "Normal Epithelial"] <- "Normal epithelial"
    lv_cts[lv_cts == "Cancer Epithelial"] <- "Cancer epithelial"

    cts <- as.character(cts)
    cts[cts == "Normal Epithelial"] <- "Normal epithelial"
    cts[cts == "Cancer Epithelial"] <- "Cancer epithelial"

    cts <- factor(cts, levels = lv_cts)
    return(cts)
}
data_melt$Row <- replace_cts(data_melt$Row)


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
        paste0(pp),
        paste0(pp),
        paste0("CGP ", pp),
        paste0("CM ", pp),
        paste0(pp)
    )
    return(pp1)
})

replace_closest_space_with_newline <- function(input_string) {
    str_length <- nchar(input_string)
    center_position <- ceiling(str_length / 2)
    space_positions <- gregexpr(" ", input_string)[[1]]
    if (length(space_positions) == 0) {
        return(input_string)
    }
    if (str_length <= 40) {
        return(input_string)
    }
    closest_space <- space_positions[which.min(abs(space_positions - center_position))]
    modified_string <- substring(input_string, 1, closest_space - 1)
    modified_string <- paste0(modified_string, "\n", substring(input_string, closest_space + 1, str_length))

    return(modified_string)
}


levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)

plot <- data_melt %>%
    filter(Pct > 1) %>%
    ggplot(aes(x = Column, y = Row, color = Score, size = Pct)) +
    geom_point() +
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
        legend.spacing.y = unit(2, "cm"),
        legend.margin = margin(t = 2.7, unit = "cm")
    )

ggsave(
    file = paste0(enrich_wd, "dotplot_cts_pathway_matrix_iCAESARunasg_pval.png"),
    plot = plot, width = 16, height = 6.5, units = "in", dpi = 200
)



## 3.3 spatial plot of each pathway on each section
features <- pathway_used
for (i in 1:4) {
    seu <- seuList[[i]]
    enrich_wd_i <- paste0(enrich_wd, "BC", i, "/")
    dir.create(enrich_wd_i, showWarnings = FALSE)

    seu <- AddMetaData(seu, as.data.frame(pathway_scoresList[[i]]))

    coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
    if (i %in% c(1:2)) {
        plot_height <- coord[2] / 500
        plot_width <- coord[1] / 500
    } else {
        plot_height <- coord[2] / 1000
        plot_width <- coord[1] / 1000
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
                legend.position = "right",
                legend.justification = "center",
                legend.box = "vertical"
            )
        ggsave(
            file = paste0(enrich_wd_i, "FeaturePlot_pos_", feature, "_BC", i, ".png"),
            plot = plot, width = plot_width, height = plot_height, units = "in", dpi = 200
        )
    }
}
