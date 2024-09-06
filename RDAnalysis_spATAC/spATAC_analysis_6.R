rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)
ana_wd <- paste0(work_wd, "spATACME11/")
setwd(ana_wd)

load("seu.rda")
seu <- pdistance(seu, reduction = "caesar", assay.name = "distce")


load("pathwayList.rda")
pathway_names_list <- lapply(pathwayList, names)
pathway_names <- Reduce(c, pathway_names_list)

load(paste0(ana_wd, "pathway_scores.rda"))


# ------------------------------------------------------------------------------
# pathways
# ------------------------------------------------------------------------------
n_path <- 5
threshold <- 0.95

gs <- pathway_names
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
ct_enrich_pval <- CAESAR.CTDEP(y, pathway_scores)

save(ct_enrich_pval, file = "ct_enrich_pval.rda")

### adjust p values
ct_enrich_pval_adj <- apply(ct_enrich_pval, 2, p.adjust, method = "BH")

## select top 5 pathway for each celltype
pathway_candidate <- rownames(ct_enrich_pval_adj)[!grepl("^GSE", rownames(ct_enrich_pval_adj))]
ntop <- 5
pathway_used <- lapply(cts, function(ct) {
    x <- ct_enrich_pval_adj[pathway_candidate, ct]
    x <- sort(x, decreasing = FALSE)
    x <- x[x <= 0.05]
    if (length(x) == 0) {
        return(NULL)
    }
    x[seq_len(min(length(x), ntop))]
})
names(pathway_used) <- cts
pathway_used <- pathway_used[!sapply(pathway_used, is.null)]

pathway_usedList <- pathway_used
save(pathway_usedList, file = "pathway_usedList.rda")

cts_used <- setdiff(cts[1:7], "unassigned")
cts_used <- cts_used[c(1, 4, 2, 3, 5, 6)]
pathway_used <- pathway_usedList
pathway_used <- lapply(pathway_used, function(p) names(p))
pathway_used <- unlist(pathway_used[cts_used], use.names = FALSE)
save(pathway_used, file = "pathway_used.rda")

cts_pathway_scores <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_used])
}))
colnames(cts_pathway_scores) <- gsub("_", " ", colnames(cts_pathway_scores))

cts_pathway_pct <- t(sapply(cts_used, function(ct) {
    colMeans(pathway_scores[y == ct, pathway_used] > threshold)
}))
colnames(cts_pathway_pct) <- gsub("_", " ", colnames(cts_pathway_pct))


data_score <- melt(cts_pathway_scores)
colnames(data_score) <- c("Row", "Column", "Score")
data_pct <- melt(cts_pathway_pct * 100)
colnames(data_pct) <- c("Row", "Column", "Pct")

data_melt <- left_join(data_score, data_pct, by = c("Row", "Column"))
data_melt$Row <- factor(data_melt$Row, levels = rownames(cts_pathway_scores))
data_melt$Column <- factor(data_melt$Column, levels = colnames(cts_pathway_scores))

save(cts_pathway_scores, cts_pathway_pct, data_melt, file = paste0(kegg_wd, "plot_data.rda"))

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtree)

gsl_all <- lapply(pathway_set_used, function(j) {
    file_name <- switch(j,
        paste0(kegg_wd, "gene.set.list_H.rda"),
        paste0(kegg_wd, "gene.set.list_C2_CPKEGG.rda"),
        paste0(kegg_wd, "gene.set.list_C2_CPREACTOME.rda"),
        paste0(kegg_wd, "gene.set.list_C2_CGP.rda"),
        paste0(kegg_wd, "gene.set.list_C3_MIR_MIRDB.rda"),
        paste0(kegg_wd, "gene.set.list_C3_TFT_GTRD.rda"),
        paste0(kegg_wd, "gene.set.list_C5_GOBP.rda"),
        paste0(kegg_wd, "gene.set.list_C5_GOCC.rda"),
        paste0(kegg_wd, "gene.set.list_C5_GOMF.rda"),
        paste0(kegg_wd, "gene.set.list_C7_IMMUNESIGDB.rda")
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
    # 计算字符串的总长度
    str_length <- nchar(input_string)

    # 找到字符串的中心位置
    center_position <- ceiling(str_length / 2)

    # 找到所有空格的位置
    space_positions <- gregexpr(" ", input_string)[[1]]

    # 如果没有空格，返回原字符串
    if (length(space_positions) == 0) {
        return(input_string)
    }
    if (str_length <= 42) {
        return(input_string)
    }

    # 找到最靠近中心位置的空格
    closest_space <- space_positions[which.min(abs(space_positions - center_position))]

    # 将最靠近中心的空格替换为换行符
    modified_string <- substring(input_string, 1, closest_space - 1)
    modified_string <- paste0(modified_string, "\n", substring(input_string, closest_space + 1, str_length))

    return(modified_string)
}

levels(data_melt$Column) <- sapply(levels(data_melt$Column), replace_closest_space_with_newline)


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
        limits = c(0, 1)
    ) +
    # scale_x_discrete(position = "top") +
    cowplot::theme_cowplot() +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),
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
        legend.text = element_text(size = 15, face = "bold"),
        legend.spacing.y = unit(5, "cm")
    )

ggsave(
    file = paste0(kegg_wd, "dotplot_cts_pathway_matrix_CAESARunasg.png"),
    plot = plot, width = 11, height = 15, units = "in", dpi = 200
)













features <- pathway_used
kegg_wd_f <- paste0(kegg_wd, "spATAC/")
dir.create(kegg_wd_f, showWarnings = FALSE)

seu <- AddMetaData(seu, as.data.frame(pathway_scores[, pathway_used]))

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
        file = paste0(kegg_wd_f, "Step5_FeaturePlot_pos_", feature, ".png"),
        plot = plot, width = 6, height = 5, units = "in", dpi = 200
    )
}






