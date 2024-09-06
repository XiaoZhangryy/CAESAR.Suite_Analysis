rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(gsdensity)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")
load(paste0(ana_wd, "seu.rda"))

# ==============================================================================
# 1. auc - fig. 4.g
# ==============================================================================
gene.used <- VariableFeatures(seu)
el <- compute.nn.edges(coembed = Embeddings(seu, "mca"))
seu <- pdistance(seu, reduction = "fastimg", assay.name = "distce")

save(el, file = paste0(ana_wd, "el_gsdenisty_edges.rda"))

# ------------------------------------------------------------------------------
# signatrue and de genes w.r.t original labels
# ------------------------------------------------------------------------------
Idents(seu) <- seu$merged_deviAnnotate
counts_ct <- sort(table(Idents(seu)), decreasing = TRUE)
Idents(seu) <- factor(Idents(seu), levels = names(counts_ct))

deg <- FindAllMarkers(seu)

save(deg, file = paste0(ana_wd, "deg_merged_deviAnnotate.rda"))

# ------------------------------------------------------------------------------
# calculate measures based on deg
# ------------------------------------------------------------------------------
degs_top3 <- deg %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)

markerList <- split(degs_top3$gene, degs_top3$cluster)

pathway_cellid <- RunCellHGT(seu, pathways = markerList, dims = 1:q_est, minSize = 2, reduction = "mca")
pathway_gsdenisty <- run.rwr.list(el = el, gene_set_list = markerList, cells = colnames(seu))
pathway_scores <- CAESAR.enrich.score(seu, markerList)

res_auc <- auc_score_List(seu$merged_deviAnnotate, list(
    pathway_scores, t(pathway_cellid), pathway_gsdenisty
), FALSE)

apply(res_auc, 1, quantile)

save(res_auc, file = paste0(ana_wd, "enrich_measures_deg.rda"))

# ------------------------------------------------------------------------------
# plot it
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)

res_auc <- t(res_auc)
colnames(res_auc) <- c("CAESAR", "Cell-ID", "GSDensity")

df <- as.data.frame(res_auc)
df$CellType <- rownames(df)
df_long <- melt(df, id.vars = "CellType", variable.name = "Method", value.name = "Value")

p_auc <- ggplot(df_long, aes(x = Method, y = Value)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.75) +
    geom_jitter(size = 2, width = 0.25, aes(fill = CellType), shape = 21, stroke = 0.5) + #
    theme_minimal() +
    ylim(c(0.4, 1)) +
    labs(title = NULL, x = NULL, y = "AUC Score") +
    scale_fill_manual(values = cols) +
    theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # Background with border
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "none"
    )
# p_auc

ggsave(
    file = paste0(figs_wd, "auc_plot.png"),
    plot = p_auc, width = 3.2, height = 4, units = "in", dpi = 200
)
