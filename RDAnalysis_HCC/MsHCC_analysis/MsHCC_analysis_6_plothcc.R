rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/MsHCC/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")
load("ms2hsgenes.rda")

seuList <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "HCC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    seu
})

################################################################################
# 1. spatial plot
################################################################################

## do it in 2

################################################################################
# 2. co-embedding plot - fig. 3.f
################################################################################
load("sgInt.rda")
load("ct_ratio.rda")
load("seu_sc.rda")
load("markerList.rda")
load("sg_List.rda")
load("sg_sc_List.rda")

cts_used <- c("HCC cell", "Macrophages", "Fibroblasts")

msmarker <- setNames(lapply(cts_used, function(ct) {
    genes <- markerList[[1]][[ct]]
    id <- sapply(genes, function(msgene) {
        hmgene <- unique(ms2hsgenes$HGNC.symbol[ms2hsgenes$MGI.symbol %in% msgene])
        sapply(hmgene, function(gene) {
            match(gene, sgInt[[ct]])
        }) %>% min()
    })
    # id[order(id, decreasing = FALSE)][1:2]
    genes[order(id, decreasing = FALSE)][1:2]
}), cts_used)

hmmarker <- lapply(msmarker, function(msgene) {
    unique(ms2hsgenes$HGNC.symbol[ms2hsgenes$MGI.symbol %in% msgene])
})

for (i in 1:4) {
    seu <- seuList[[i]]
    Idents(seu) <- seu$CAESARunasg
    ana_wd <- paste0(work_wd, "HCC", i, "/")

    ct_used <- names(ct_ratio[[i]][ct_ratio[[i]] > 0.075])
    ct_used <- intersect(ct_used, celltypes)

    df_coumap_gene <- Reduce(rbind, lapply(
        ct_used, function(ct) {
            sig_sp <- sg_List[[i]][[ct]]
            sig_sp[match(hmmarker[[ct]], sig_sp$gene), ]
        }
    ))
    coumap_gene <- df_coumap_gene %>%
        .$gene %>%
        unique()

    seu <- CoUMAP(
        seu,
        reduction = "caesar", reduction.name = "CoUMAP",
        gene.set = coumap_gene, n_neighbors = 20, min_dist = 0.1, spread = 0.4
    )

    plot3 <- CoUMAP.plot(
        seu,
        reduction = "CoUMAP", gene_txtdata = df_coumap_gene,
        cols = c("gene" = "#000000", cols),
        pt_size = 0.6, pt_text_size = 4, alpha = 0.8
    ) + guides(
        shape = guide_legend(override.aes = list(size = 3)),
        color = guide_legend(ncol = 1, override.aes = list(size = 3))
    )
    ggsave(
        file = paste0(ana_wd, "CoUMAP_CAESARunasg_HCC", i, "_marker.png"),
        plot = plot3, width = 9, height = 5, units = "in", dpi = 200
    )
}



df_coumap_gene_sc <- Reduce(rbind, lapply(
    cts_used, function(ct) {
        sig_sp <- sg_sc_List[[1]][[ct]]
        sig_sp[match(msmarker[[ct]], sig_sp$gene), ]
    }
))

coumap_gene_sc <- df_coumap_gene_sc %>%
    .$gene %>%
    unique()

Idents(seu_sc) <- factor(seu_sc$merged_clusters, levels = names(cols))

seu_sc <- CoUMAP(
    seu_sc, reduction = "ncfm", reduction.name = "CoUMAP",
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
    file = paste0(figs_wd, "sc_CoUMAP_CAESARunasg_marker.png"),
    plot = plot3, width = 9, height = 5, units = "in", dpi = 200
)

