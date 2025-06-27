rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")
ana_wd <- paste0(work_wd, "spATACME/")

load("cols_cts_q.rda")

setwd(ana_wd)

# ------------------------------------------------------------------------------
# pathways
# ------------------------------------------------------------------------------
n_path <- 5
threshold <- 0.95

gsl <- lapply(pathway_set_used, function(j) {
    file_name <- switch(j,
        paste0(enrich_wd, "gene.set.list_H.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CPKEGG.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CPREACTOME.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CGP.rda"),
        paste0(enrich_wd, "gene.set.list_C3_MIR_MIRDB.rda"),
        paste0(enrich_wd, "gene.set.list_C3_TFT_GTRD.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOBP.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOCC.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOMF.rda"),
        paste0(enrich_wd, "gene.set.list_C7_IMMUNESIGDB.rda")
    )
    load(file_name)
    names(gene.set.list)
})
gs <- unlist(gsl)
gs <- gs[!grepl("^GSE", gs)]

y <- as.character(seu$CAESARunasg)
ct_count <- sort(table(y), decreasing = TRUE)
cts <- names(ct_count)
cps <- sapply(cts, function(ct) {
    colMeans(pathway_scores[y == ct, ])
})
cps_pt <- sapply(cts, function(ct) {
    colMeans(pathway_scores[y == ct, ] > threshold)
})

## calculate celltype enrichment
ct_enrich <- function(y, pathway_scores, cts = NULL) {
    if (is.null(cts)) {
        cts <- names(sort(table(y), decreasing = TRUE))
    }
    sapply(cts, function(ct) {
        id1 <- which(y == ct)
        id2 <- which(y != ct)
        apply(pathway_scores, 2, function(x) {
            wilcox.test(x[id1], x[id2], alternative = "greater")$p.value
        })
    })
}
ct_enrich_pval <- ct_enrich(y, pathway_scores)

save(ct_enrich_pval, file = "ct_enrich_pval.rda")

### adjust p values
ct_enrich_pval_adj <- apply(ct_enrich_pval, 2, p.adjust, method = "BH")

## select top 5 pathway for each celltype
pathway_candidate <- rownames(ct_enrich_pval_adj)[!grepl("^GSE", rownames(ct_enrich_pval_adj))]
ntop <- 5
pathway_usedList <- lapply(cts, function(ct) {
    x <- ct_enrich_pval_adj[pathway_candidate, ct]
    x <- sort(x, decreasing = FALSE)
    x <- x[x <= 0.05]
    if (length(x) == 0) {
        return(NULL)
    }
    x[seq_len(min(length(x), ntop))]
})
names(pathway_usedList) <- cts
pathway_usedList <- pathway_usedList[!sapply(pathway_usedList, is.null)]

# save(pathway_usedList, file = "pathway_usedList.rda")

cts_used <- setdiff(cts[1:7], "unassigned")
cts_used <- cts_used[c(1, 4, 2, 3, 5, 6)]
pathway_usedList <- lapply(pathway_usedList, function(p) names(p))
pathway_usedList <- unlist(pathway_usedList[cts_used], use.names = FALSE)
save(pathway_usedList, file = "pathway_usedList.rda")










load("pathway_usedList.rda")
cts_pathway_scores <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_usedList])
}))
colnames(cts_pathway_scores) <- gsub("_", " ", colnames(cts_pathway_scores))

cts_pathway_pct <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_usedList] > threshold)
}))
colnames(cts_pathway_pct) <- gsub("_", " ", colnames(cts_pathway_pct))


data_score <- melt(cts_pathway_scores)
colnames(data_score) <- c("Row", "Column", "Score")
data_pct <- melt(cts_pathway_pct * 100)
colnames(data_pct) <- c("Row", "Column", "Pct")

data_melt <- left_join(data_score, data_pct, by = c("Row", "Column"))
data_melt$Row <- factor(data_melt$Row, levels = rownames(cts_pathway_scores))
data_melt$Column <- factor(data_melt$Column, levels = colnames(cts_pathway_scores))

save(cts_pathway_scores, cts_pathway_pct, data_melt, file = paste0(enrich_wd, "plot_data.rda"))

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtree)

gsl_all <- lapply(pathway_set_used, function(j) {
    file_name <- switch(j,
        paste0(enrich_wd, "gene.set.list_H.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CPKEGG.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CPREACTOME.rda"),
        paste0(enrich_wd, "gene.set.list_C2_CGP.rda"),
        paste0(enrich_wd, "gene.set.list_C3_MIR_MIRDB.rda"),
        paste0(enrich_wd, "gene.set.list_C3_TFT_GTRD.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOBP.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOCC.rda"),
        paste0(enrich_wd, "gene.set.list_C5_GOMF.rda"),
        paste0(enrich_wd, "gene.set.list_C7_IMMUNESIGDB.rda")
    )
    load(file_name)
    names(gene.set.list)
})

levels(data_melt$Column) <- sapply(colnames(cts_pathway_scores), function(pp) {
    pp0 <- gsub(" ", "_", pp)
    j <- which(sapply(seq_along(gsl_all), function(j) {
        pp0 %in% gsl_all[[j]]
    }))
    pp <- gsub("MODULE", "CANCER RELATED MODULE", pp)
    pp1 <- switch(pathway_set_used[j],
        paste0("H ", pp),
        paste0("KEGG ", pp),
        paste0(pp),
        paste0("CGP ", pp),
        paste0("MIRDB ", pp),
        paste0("GTRD ", pp),
        paste0(pp),
        paste0(pp),
        paste0(pp),
        paste0("IMMUNESIGDB", pp)
    )
    return(pp1)
})

max(sapply(levels(data_melt$Column), nchar) / 2)

replace_closest_space_with_newline <- function(input_string) {
    str_length <- nchar(input_string)
    center_position <- ceiling(str_length / 2)
    space_positions <- gregexpr(" ", input_string)[[1]]
    if (length(space_positions) == 0) {
        return(input_string)
    }
    if (str_length <= 42) {
        return(input_string)
    }
    closest_space <- space_positions[which.min(abs(space_positions - center_position))]
    modified_string <- substring(input_string, 1, closest_space - 1)
    modified_string <- paste0(modified_string, "\n", substring(input_string, closest_space + 1, str_length))

    return(modified_string)
}

levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)



### Fig 5.d
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
        limits = c(0, 1)
    ) +
    cowplot::theme_cowplot() +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 14),
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
        legend.spacing.y = unit(3, "cm")
    )

ggsave(
    file = paste0(enrich_wd, "dotplot_cts_pathway_matrix_CAESARunasg.png"),
    plot = plot, width = 20, height = 9, units = "in", dpi = 200
)










### Fig 5.e


features <- pathway_usedList
enrich_wd_f <- paste0(enrich_wd, "pathways/")
dir.create(enrich_wd_f, showWarnings = FALSE)

seu <- AddMetaData(seu, as.data.frame(pathway_scores[, pathway_usedList]))

for (feature in features) {
    expr_values <- FetchData(seu, vars = feature)

    plot <- FeaturePlot(seu, features = feature, reduction = "pos") +
        scale_color_gradientn(
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
        file = paste0(enrich_wd_f, "Step5_FeaturePlot_pos_", feature, ".png"),
        plot = plot, width = 6, height = 5, units = "in", dpi = 200
    )
}














