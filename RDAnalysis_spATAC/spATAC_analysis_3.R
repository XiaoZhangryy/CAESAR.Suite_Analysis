rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)


################################################################################
# 1. load data
################################################################################
## 1.1 scRNA
load("/share/analysisdata/zhangx/CAESAR/SpATACME/seu_sc_raw.rda")

## 1.2 spatac
load("/share/analysisdata/zhangx/CAESAR/SpATACME/spATAC/seuATAC_final.rda")

seu_sp <- seuATAC
rm(seuATAC)

print(seu_sc)
print(seu_sp)

################################################################################
# 1. preprocess
################################################################################
work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

cols <- setNames(
    c(
        "#034e7b", "#3690c0", "#a6bddb", "#7AC3DF", "#0570b0", # cts_1
        "#70CDBE", "#00441b", "#006d2c", "#41ab5d", "#74c476", "#a1d99b",
        "#c7e9c0", "#e5f5e0", "#238b45", # cts_2
        "#addd8e", "#67000d", "#cb181d", "#fc9272", "#fcbba1",
        "#005a32", "#dd3497", "#fa9fb5", "#fcc5c0", "#fde0dd",
        "#a63603", "#fec44f", "#238443", "#fee391", "#ffff99",
        "#737373"
    ),
    c(
        "Chondrocytes & osteoblasts", "Intermediate Mesoderm",
        "Chondroctye progenitors", "Connective tissue progenitors",
        "Osteoblasts", # cts_1
        "Excitatory neurons", "Postmitotic premature neurons",
        "Isthmic organizer cells", "Inhibitory interneurons",
        "Inhibitory neuron progenitors", "Inhibitory neurons",
        "Neural progenitor cells", "Oligodendrocyte Progenitors",
        "Radial glia", # cts_2
        "Neural Tube", "Epithelial cells", "Limb mesenchyme",
        "Premature oligodendrocyte", "Early mesenchyme", "Sensory neurons",
        "Endothelial cells",
        "Stromal cells", "Notochord cells", "Primitive erythroid lineage",
        "Definitive erythroid lineage", "Ependymal cell",
        "Granule neurons", "White blood cells", "Hepatocytes", # cts_3
        "unassigned"
    )
)
celltypes <- setdiff(names(cols), "unassigned")
q_est <- 50
save(cols, celltypes, q_est, file = "cols_cts_q.rda")

## -----------------------------------------------------------------------------
## 1.1 align scRNA and spatac data
## -----------------------------------------------------------------------------
common_genes <- intersect(rownames(seu_sc), rownames(seu_sp))
print(length(common_genes))
# [1] 18155

## 1.1.1 align genes
seu_sc <- seu_sc[common_genes, ]
seu_sp <- seu_sp[common_genes, ]

### 1.1.2 Normalize
seu_sc <- seu_sc %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000)
seu_sp <- seu_sp %>%
    FindVariableFeatures(nfeatures = 2000)

### 1.1.3 align variable genes
common_vf <- VariableFeatures(seu_sc)

VariableFeatures(seu_sc) <- common_vf
VariableFeatures(seu_sp) <- common_vf


################################################################################
# 2. extract signature genes on reference
################################################################################
## 2.1 co-embedding
seu_sc <- NCFM(seu_sc, q = q_est)

save(seu_sc, file = "seu_sc.rda")

## 2.2 find sg
seu_sc <- pdistance(seu_sc, reduction = "ncfm")
sg_sc_List <- find.sig.genes(seu_sc) %>% list()
save(sg_sc_List, file = "sg_sc_List.rda")

### 2.2.5 find marker
markerList <- lapply(sg_sc_List, function(sig) {
    marker <- marker.select(sig, overlap.max = 3)
    marker
})

n.marker <- sapply(markerList, function(marker) length(marker[[1]]))
message("The number of markers for each reference dataset are (", paste0(n.marker, collapse = ", "), ")")

save(markerList, file = "markerList.rda")

################################################################################
# 3. save spatial data
################################################################################
cords <- c("pxl_row_in_fullres", "pxl_col_in_fullres")
pos <- as.matrix(seu_sp@meta.data[, cords])
pos[, 2] <- -pos[, 2]
colnames(pos) <- paste0("pos", 1:2)

seu_sp@reductions[["pos"]] <- CreateDimReducObject(
    embeddings = pos,
    key = paste0("pos", "_"), assay = "RNA"
)

seu <- seu_sp
save(seu, file = "seu_processed.rda")

print(apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min))

print(seu_sc)
print(seu)


################################################################################
# 4. enrichment
################################################################################
library(msigdbr)

gene.used <- VariableFeatures(seu_sc)

pathwayList <- lapply(1:10, function(j) {
    ms_KEGG <- switch(j,
        msigdbr(species = "Mus musculus", category = "H", subcategory = ""),
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG"),
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP"),
        msigdbr(species = "Mus musculus", category = "C3", subcategory = "MIR:MIRDB"),
        msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF"),
        msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")
    )

    pathway_list <- ms_KEGG %>%
        group_by(gs_name) %>%
        summarise(genes = list(intersect(gene_symbol, gene.used))) %>%
        tibble::deframe()
    n.pathway_list <- sapply(pathway_list, length)
    pathway_list <- pathway_list[n.pathway_list >= 5]

    pathway_list
})
save(pathwayList, file = "pathwayList.rda")







plot1 <- DimPlot(seu_sc, reduction = "umap_ncfm", cols = cols, label = TRUE) +
    theme(legend.position = "right") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))
ggsave(
    file = paste0(figs_wd, "sc_DimPlot_", "umap_ncfm", "_merged_clusters.png"),
    plot = plot1, width = 15, height = 10, units = "in", dpi = 200
)

plot1 <- DimPlot(seu_sc, reduction = "umap_mca", cols = cols, label = TRUE) +
    theme(legend.position = "right") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))
ggsave(
    file = paste0(figs_wd, "sc_DimPlot_", "umap_mca", "_merged_clusters.png"),
    plot = plot1, width = 15, height = 10, units = "in", dpi = 200
)

plot1 <- DimPlot(seu, reduction = "pos", group.by = "Clusters_unsupervised") +
    theme(
        axis.line = element_blank(), # Remove the axis lines
        axis.text.x = element_blank(), # Remove the text on the x-axis
        axis.text.y = element_blank(), # Remove the text on the y-axis
        axis.ticks = element_blank(), # Remove the ticks
        axis.title.x = element_blank(), # Remove the x-axis title
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold")
    ) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))
ggsave(
    file = paste0(figs_wd, "sp_DimPlot_", "pos", "_Clusters_unsupervised.png"),
    plot = plot1, width = 10, height = 5, units = "in", dpi = 200
)


plot1 <- DimPlot(
    seu,
    reduction = "pos", group.by = "Clusters_unsupervised",
    split.by = "Clusters_unsupervised", ncol = 2
) +
    theme(
        axis.line = element_blank(), # Remove the axis lines
        axis.text.x = element_blank(), # Remove the text on the x-axis
        axis.text.y = element_blank(), # Remove the text on the y-axis
        axis.ticks = element_blank(), # Remove the ticks
        axis.title.x = element_blank(), # Remove the x-axis title
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold")
    ) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))
ggsave(
    file = paste0(figs_wd, "sp_DimPlot_", "pos", "_split_Clusters_unsupervised.png"),
    plot = plot1, width = 21, height = 12, units = "in", dpi = 200
)

seu <- RunUMAP(
    seu,
    reduction = "lsi", dims = 1:30, reduction.name = "umap_lsi"
)

for (redu in c("umap_lsi", "pos")) {
    for (feat in c("predictedGroup")) {
        plot1 <- DimPlot(seu, reduction = redu, group.by = feat, cols = cols) +
            theme(
                axis.line = element_blank(), # Remove the axis lines
                axis.text.x = element_blank(), # Remove the text on the x-axis
                axis.text.y = element_blank(), # Remove the text on the y-axis
                axis.ticks = element_blank(), # Remove the ticks
                axis.title.x = element_blank(), # Remove the x-axis title
                axis.title.y = element_blank(),
                legend.position = "right",
                legend.text = element_text(size = 10, face = "bold")
            ) +
            guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))
        ggsave(
            file = paste0(figs_wd, "sp_DimPlot_", redu, "_", feat, ".png"),
            plot = plot1, width = 10, height = 5, units = "in", dpi = 200
        )


        plot1 <- DimPlot(
            seu,
            reduction = redu, group.by = feat,
            split.by = feat, ncol = 5, cols = cols
        ) +
            theme(
                axis.line = element_blank(), # Remove the axis lines
                axis.text.x = element_blank(), # Remove the text on the x-axis
                axis.text.y = element_blank(), # Remove the text on the y-axis
                axis.ticks = element_blank(), # Remove the ticks
                axis.title.x = element_blank(), # Remove the x-axis title
                axis.title.y = element_blank(),
                legend.position = "right",
                legend.text = element_text(size = 10, face = "bold")
            ) +
            guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))
        ggsave(
            file = paste0(figs_wd, "sp_DimPlot_", redu, "_split_", feat, ".png"),
            plot = plot1, width = 21, height = 12, units = "in", dpi = 200
        )
    }
}





# cell embedding
library(Seurat)
library(data.table)
library(ProFAST)
library(dplyr)
library(ggplot2)
library(Seurat)
library(CelliD)
library(DR.SC)

dir_image_feature <- "/share/analysisdata/liuw/coembed_image/RealData/MouseEmbryo/ME11/"
library(CoFAST)

meta_data <- read.csv("/share/analysisdata/liuw/coembed_image/RealData/MouseEmbryo/ME11/meta_data.csv")
## Check the bardcode aligned
print(all(meta_data$barcode == colnames(seu)))

set.seed(2024)
feature_img <- fread(paste0(dir_image_feature, "feature_img.csv"))
feature_img <- as.matrix(feature_img)

prin <- DR.SC:::wpca(feature_img, q = 10, weighted = FALSE)
embed_img <- prin$PCs
metadata <- seu@meta.data
pos <- as.matrix(metadata[c("array_row", "array_column")])
## Here, the neighbors are four

radius_use <- searchRadius(pos, lower.med = 3.5, upper.med = 5.5, radius.upper = 1.8)


set.seed(1)
n_spots <- nrow(pos)
idx <- sample(n_spots, min(100, n_spots))
dis <- dist(embed_img[idx, ])
sigma <- median(dis)^2
wAdj <- weightAdj(pos, img_embed = embed_img, radius = radius_use, width = sigma)
Matrix::colSums(wAdj)
Matrix::colSums(wAdj > 0)
wAdj[1, wAdj[1, ] > 0]

colnames(seu)[1:4]
seu <- FAST_IMG(seu, weightAdj = wAdj, q = q_est, approx_Phi = TRUE)





seu_name <- paste0(ana_wd, "seu.rda")
save(seu, file = seu_name)

seu_sc@assays[["distce"]] <- NULL
seu_sc_name <- paste0(ana_wd, "seu_sc.rda")
save(seu_sc, file = seu_sc_name)
