rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/PixelMOB/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

cols <- setNames(
    c(
        "#7A48A4", "#5DADE2", "#d9f0a3", "#00441b", "#ffffb2",
        "#D8FDF7", "#fee391", "#E04D50", "#EA6F69", "#1A5276",
        "#f768a1", "#FBCCE2", "#aa8b00", "#FFA500", "#807dba",
        "#c7e9c0", "#bcbddc", "#E3BF9F", "#737373"
    ),
    c(
        "Astro", "EC", "ImmunoCells", "Mes", "Microglia",
        "Mural", "MyOligo", "Neuron.Astro-Like", "EPL-IN", "GC",
        "Neuron.Immature", "M/TC", "Neuron.OSN", "PGC", "Neuron.Transition",
        "OEC", "OPC", "RBCs", "unassigned"
    )
)
celltypes <- setdiff(names(cols), "unassigned")
q_est <- 50
save(cols, celltypes, q_est, file = "cols_cts_q.rda")

################################################################################
# 1. preprocess
################################################################################
## -----------------------------------------------------------------------------
## 1.1 read scRNA data
## -----------------------------------------------------------------------------
read_scRNA_PixelMOB <- function() {
    data_path <- "/share/rawdata/zhangx/scRNA/mouse_OB/GSE121891/"
    setwd(data_path)

    OB_metaData_seurat <- read.csv("GSE121891_OB_metaData_seurat.csv")
    OB_6_runs.raw.dge <- read.csv("GSE121891_OB_6_runs.raw.dge.csv")
    OB_6_runs_processed_seurat.dge <- read.csv("GSE121891_OB_6_runs_processed_seurat.dge.csv")
    Figure_2_metadata <- read.table("GSE121891_Figure_2_metadata.txt", header = TRUE, sep = "", dec = ".")
    rownames(OB_metaData_seurat) <- OB_metaData_seurat$X
    OB_metaData_seurat$X <- NULL

    FinalCelltype <- OB_metaData_seurat$ClusterName
    names(FinalCelltype) <- rownames(OB_metaData_seurat)
    FinalCelltype[rownames(Figure_2_metadata)] <- Figure_2_metadata$`FinalIds`

    OB_metaData_seurat$`res.1` <- NA
    OB_metaData_seurat$`ClusterId` <- NA
    OB_metaData_seurat$`stim` <- NA
    OB_metaData_seurat$`celltype.stim` <- NA
    OB_metaData_seurat$`FinalIds` <- NA
    OB_metaData_seurat[rownames(Figure_2_metadata), "res.1"] <- Figure_2_metadata$`res.1`
    OB_metaData_seurat[rownames(Figure_2_metadata), "ClusterId"] <- Figure_2_metadata$`ClusterId`
    OB_metaData_seurat[rownames(Figure_2_metadata), "stim"] <- Figure_2_metadata$`stim`
    OB_metaData_seurat[rownames(Figure_2_metadata), "celltype.stim"] <- Figure_2_metadata$`celltype.stim`
    OB_metaData_seurat[rownames(Figure_2_metadata), "FinalIds"] <- Figure_2_metadata$`FinalIds`
    OB_metaData_seurat$FinalCelltype <- FinalCelltype

    print(dim(OB_6_runs_processed_seurat.dge))

    seu_sc <- CreateSeuratObject(
        counts = OB_6_runs_processed_seurat.dge,
        meta.data = OB_metaData_seurat, min.cells = 1, min.features = 5
    )

    print(dim(seu_sc))

    seu_sc <- subset(seu_sc, subset = FinalCelltype %in% c("N1", "N2", "N3", "N8", "N9", "N11", "N14"), invert = TRUE)

    print(dim(seu_sc))

    Idents(seu_sc) <- seu_sc$FinalCelltype

    newcluster.ids <- c(
        "Astro", "Astro", "Astro", "EC", "EC", "Mes", "Mes",
        "Microglia", "Microglia", "Microglia", "ImmunoCells",
        "Mural", "Mural", "MyOligo", "ImmunoCells", "Neuron.OSN", "PGC",
        "GC", "Neuron.Immature", "PGC", "Neuron.Transition",
        "GC", "PGC", "GC", "GC", "GC",
        "GC", "Neuron.Astro-Like", "GC", "M/TC",
        "M/TC", "M/TC", "EPL-IN", "OEC", "OEC", "OEC",
        "OEC", "OEC", "OPC", "RBCs"
    )
    names(newcluster.ids) <- sort(levels(seu_sc))
    seu_sc <- RenameIdents(seu_sc, newcluster.ids)
    seu_sc$merged_clusters <- Idents(seu_sc)

    print(dim(seu_sc))

    print(newcluster.ids)

    setwd(work_wd)

    return(seu_sc)
}

seu_sc <- read_scRNA_PixelMOB()

table(seu_sc$merged_clusters)

## -----------------------------------------------------------------------------
## 1.2 read Xenium data
## -----------------------------------------------------------------------------
read_sp_PixelMOB <- function() {
    data_path <- "/share/rawdata/zhangx/Pixel_seq/mouse_OB_PBN/GSE186097/"
    seu_sp <- readRDS(paste0(data_path, "GSM5631821_OB_seurat.rds"))

    OB36_metadata <- read.csv(paste0(data_path, "OB36_metadata.csv"))
    rownames(OB36_metadata) <- OB36_metadata$X
    OB36_metadata$X <- NULL

    pos <- as.matrix(seu_sp@images$image@coordinates[, c("x", "y")])
    colnames(pos) <- paste0("pos_", 1:2)
    seu_sp@reductions[["pos"]] <- CreateDimReducObject(
        embeddings = pos, key = "pos_", assay = "RNA"
    )

    seu_sp <- AddMetaData(
        object = seu_sp,
        metadata = OB36_metadata[, c("deviAnnotate", "annotateName")],
        col.name = c("deviAnnotate", "annotateName")
    )

    seu_sp$x <- Embeddings(seu_sp, "pos")[, 1]
    seu_sp$y <- Embeddings(seu_sp, "pos")[, 2]

    print(dim(seu_sp[["Spatial"]]))

    seu <- CreateSeuratObject(
        seu_sp[["Spatial"]]@counts,
        project = "OB_Pixelseq",
        meta.data = seu_sp@meta.data, min.features = 5, min.cells = 1
    )
    Idents(seu) <- seu$deviAnnotate

    cords <- c("x", "y")
    pos <- as.matrix(seu@meta.data[, cords])
    colnames(pos) <- paste0("pos", 1:2)

    seu@reductions[["pos"]] <- CreateDimReducObject(
        embeddings = pos,
        key = paste0("pos", "_"), assay = "RNA"
    )

    newcluster.ids <- c(
        "Astro", "EC", "ImmunoCells", "Mes", "Microglia",
        "Mural", "MyOligo", "PGC", "GC", "PGC", "Neuron.Transition", "GC",
        "PGC", "GC", "GC", "GC", "GC", "GC", "EPL-IN", "Neuron.Astro-Like",
        "Neuron.Immature", "M/TC", "Neuron.OSN", "OEC", "OPC", "RBCs"
    )
    names(newcluster.ids) <- sort(levels(seu))
    seu <- RenameIdents(seu, newcluster.ids)
    seu$merged_deviAnnotate <- Idents(seu)

    print(dim(seu))
    print(newcluster.ids)

    seu
}

seu <- read_sp_PixelMOB()

## -----------------------------------------------------------------------------
## 1.3 align scRNA and Xenium data
## -----------------------------------------------------------------------------
common_genes <- intersect(rownames(seu_sc), rownames(seu))
print(length(common_genes))

### 1.3.1 align genes
seu_sc <- seu_sc[common_genes, ]
seu <- seu[common_genes, ]

### 1.3.2 Normalize
seu_sc <- seu_sc %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3500)
seu <- seu %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3500)

### 1.3.3 align variable genes
common_vgs <- intersect(VariableFeatures(seu_sc), VariableFeatures(seu))
print(length(common_vgs))

VariableFeatures(seu_sc) <- common_vgs
VariableFeatures(seu) <- common_vgs


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

## 2.3 find marker
markerList <- lapply(sg_sc_List, function(sig) {
    marker <- marker.select(sig, overlap.max = 2)
    marker
})

n.marker <- sapply(markerList, function(marker) length(marker[[1]]))
message("The number of markers for each reference dataset are (", paste0(n.marker, collapse = ", "), ")")

save(markerList, file = "markerList.rda")

################################################################################
# 3. save spatial data
################################################################################
save(seu, file = "seu_processed.rda")

apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)

################################################################################
# 4. enrichment
################################################################################
library(msigdbr)

gene.used <- VariableFeatures(seu_sc)

pathwayList <- lapply(1:7, function(j) {
    hm_KEGG <- switch(j,
        msigdbr(species = "Mus musculus", category = "H", subcategory = ""),
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG"),
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
        msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC"),
        msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF")
    )

    pathway_list <- hm_KEGG %>%
        group_by(gs_name) %>%
        summarise(genes = list(intersect(gene_symbol, gene.used))) %>%
        tibble::deframe()
    n.pathway_list <- sapply(pathway_list, length)
    pathway_list <- pathway_list[n.pathway_list >= 5]

    pathway_list
})
save(pathwayList, file = "pathwayList.rda")

################################################################################
# 5. Plot UMI density - fig. 4.c left panel
################################################################################
# seu$log2UMI <- log2(seu$nCount_RNA)
seu$logUMI <- log(seu$nCount_RNA)

# Create a custom scale function
scale_color_custom <- function(...) {
    scale_color_gradientn(
        colours = c("#f6eff7", "#feebe2", "#f768a1", "#7a0177", "#6e016b"),
        values = scales::rescale(c(5, 7, 8, 9, 10)),
        limits = c(min(seu$logUMI), max(seu$logUMI)), ...
    )
}

coord <- apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)
plot_height <- coord[2] / 2000
plot_width <- coord[1] / 2000

UMI_density_plot <- FeaturePlot(
    seu,
    features = "logUMI", reduction = "pos", pt.size = 0.4
) +
    NoLegend() +
    theme(
        axis.line = element_blank(), # Remove the axis lines
        axis.text.x = element_blank(), # Remove the text on the x-axis
        axis.text.y = element_blank(), # Remove the text on the y-axis
        axis.ticks = element_blank(), # Remove the ticks
        axis.title.x = element_blank(), # Remove the x-axis title
        axis.title.y = element_blank(),
        plot.title = element_blank()
    ) +
    scale_color_custom()
ggsave(
    file = paste0(figs_wd, "UMI_density.png"),
    plot = UMI_density_plot, width = plot_width, height = plot_height, units = "in", dpi = 500
)
