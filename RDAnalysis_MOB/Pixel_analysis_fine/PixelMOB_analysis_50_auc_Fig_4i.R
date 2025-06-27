rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "PixelMOB/")
load(paste0(ana_wd, "seu.rda"))

# ==============================================================================
# 1. auc
# ==============================================================================
library(gsdensity)
el <- compute.mca(object = seu, dims.use = 1:q_est) %>% compute.nn.edges()

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
# markerList <- find_marker(list(deg), 5, 1)

# save(markerList, file = "auc_marker.rda")

load(paste0(ana_wd, "auc_marker.rda"))

pathway_cellid <- RunCellHGT(seu, pathways = markerList, dims = 1:q_est, minSize = 2, reduction = "mca")
pathway_gsdenisty <- run.rwr.list(el = el, gene_set_list = markerList, cells = colnames(seu))

pathway_scores <- CAESAR.enrich.score(seu, pathwaylist = markerList)

pathway_gsea <- GSEA_enrich(seu, pathways = markerList)

pathway_gsva <- GSVA_enrich(seu, pathways = markerList)

pathway_aucell <- AUCell_enrich(seu, pathways = markerList)

pathway_vam <- VAM_enrich(seu, pathways = markerList)


res_auc <- Reduce(cbind, list(
    "CAESAR" = calculate_auc(seu$merged_deviAnnotate, pathway_scores, return.mean = FALSE),
    "Cell-ID" = calculate_auc(seu$merged_deviAnnotate, as.matrix(t(pathway_cellid)), return.mean = FALSE),
    "GSDensity" = calculate_auc(seu$merged_deviAnnotate, as.matrix(pathway_gsdenisty), return.mean = FALSE),
    "GSEA" = calculate_auc(seu$merged_deviAnnotate, pathway_gsea, return.mean = FALSE),
    "GSVA" = calculate_auc(seu$merged_deviAnnotate, pathway_gsva, return.mean = FALSE),
    "AUCell" = calculate_auc(seu$merged_deviAnnotate, pathway_aucell, return.mean = FALSE),
    "VAM" = calculate_auc(seu$merged_deviAnnotate, pathway_vam, return.mean = FALSE)
))

save(res_auc, file = paste0(ana_wd, "enrich_measures_deg.rda"))


# ------------------------------------------------------------------------------
# plot it
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)

df <- as.data.frame(res_auc)
df$CellType <- rownames(df)
df_long <- melt(df, id.vars = "CellType", variable.name = "Method", value.name = "Value")
df_long$Method <- factor(df_long$Method, levels = c(
  "CAESAR", "AUCell", "VAM", "GSVA", "GSEA", "GSDensity", "Cell-ID"
))

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
    plot = p_auc, width = 5, height = 4, units = "in", dpi = 200
)
