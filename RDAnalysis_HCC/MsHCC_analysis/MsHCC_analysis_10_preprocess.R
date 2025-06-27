rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/MsHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

cols <- setNames(
    c(
        "#fdc086", "#7fc97f", "#ffff99", "#b30000", "#a6cee3",
        "#8c6bb1", "#386cb0", "#737373"
    ),
    c(
        "B/Plasma", "Cholangiocytes", "Endothelial", "HCC cell",
        "Macrophages", "Neutrophil", "Fibroblasts", "unassigned"
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
read_scRNA_msHCC <- function() {
    seu_sc <- readRDS("/share/rawdata/zhangx/scRNA/mouse_HCC/6_samples_integrated_f500_npcs50_res2.rds")

    print(dim(seu_sc))

    seu_sc$Cluster_type <- paste(seu_sc$type, seu_sc$seurat_clusters, sep = "_")
    Idents(seu_sc) <- seu_sc$Cluster_type

    # names from https://github.com/holab-hku/HCC_Stem_scRNAseq/blob/main/scripts/Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation/README_1_standard_workflow_before_rename.R line 328
    newcluster.ids <- c("Macrophages/Monocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Endothelial_cells", "HCC/Hepatocytes", "HCC/Hepatocytes", "Endothelial_cells", "Cholangiocytes", "HCC/Hepatocytes", "Endothelial_cells", "HCC/Hepatocytes", "Cholangiocytes", "HCC/Hepatocytes", "Macrophages/Monocytes", "Macrophages/Monocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Endothelial_cells", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Stellate_cells/Fibroblasts", "B_cells/Plasma_cells", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Neutrophil", "HCC/Hepatocytes", "Neutrophil", "Endothelial_cells", "Macrophages/Monocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Cholangiocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "HCC/Hepatocytes", "Stellate_cells/Fibroblasts", "Endothelial_cells", "Endothelial_cells", "Endothelial_cells", "Endothelial_cells", "Endothelial_cells", "Endothelial_cells", "B_cells/Plasma_cells", "Endothelial_cells")
    names(newcluster.ids) <- levels(seu_sc)
    seu_sc <- RenameIdents(seu_sc, newcluster.ids)
    seu_sc$merged_clusters1 <- Idents(seu_sc)

    lineage.ids <- c(
        "Macrophages", "HCC cell", "Endothelial", "Cholangiocytes",
        "Fibroblasts", "B/Plasma", "Neutrophil"
    )
    names(lineage.ids) <- levels(seu_sc)
    seu_sc <- RenameIdents(seu_sc, lineage.ids)
    seu_sc$merged_clusters <- Idents(seu_sc)

    seu_sc <- subset(x = seu_sc, subset = day_type %in% c("Day3_Neg", "Day10_Neg", "Day30_Neg"))

    print(dim(seu_sc))

    seu_sc <- CreateSeuratObject(
        counts = seu_sc@assays$RNA@counts, meta.data = seu_sc@meta.data,
        min.features = 5, min.cells = 1
    )

    print(dim(seu_sc))

    Idents(seu_sc) <- seu_sc$merged_clusters

    return(seu_sc)
}

seu_sc <- read_scRNA_msHCC()

table(seu_sc$merged_clusters)

### mouse to human

library(biomaRt)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

msgenes <- rownames(seu_sc)
ms2hsgenes <- getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = msgenes,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = TRUE
)
save(ms2hsgenes, file = "ms2hsgenes.rda")


## -----------------------------------------------------------------------------
## 1.2 read Xenium data
## -----------------------------------------------------------------------------
read_sp_msHCC <- function() {
    # HCC data copy from /share/home/chenjy/SDcon/HCC/HCC1_seu.RDS

    seuList <- setNames(lapply(1:4, function(i) {
        readRDS(paste0("/share/analysisdata/zhangx/HCC_CoFAST/rawData/HCC", i, "_seu.RDS"))
    }), paste0("HCC", 1:4))

    print(sapply(seuList, dim))

    seuList <- lapply(seuList, function(seu) {
        seu <- CreateSeuratObject(
            counts = seu@assays$RNA@counts, meta.data = seu@meta.data,
            min.features = 5, min.cells = 1
        )

        cords <- c("imagecol", "imagerow")
        pos <- as.matrix(seu@meta.data[, cords])
        pos[, 2] <- 3500 - pos[, 2]
        colnames(pos) <- paste0("pos", 1:2)

        seu@reductions[["pos"]] <- CreateDimReducObject(
            embeddings = pos,
            key = paste0("pos", "_"), assay = "RNA"
        )
        seu
    })

    print(sapply(seuList, dim))

    seuList
}

seuList <- read_sp_msHCC()

## -----------------------------------------------------------------------------
## 1.3 align scRNA and hcc data
## -----------------------------------------------------------------------------
load("ms2hsgenes.rda")

common_genes <- Reduce(intersect, c(
    list(ms2hsgenes$HGNC.symbol),
    lapply(seuList, rownames)
))
print(length(common_genes))

### 1.3.1 align genes
mouse_common_genes <- unique(ms2hsgenes$MGI.symbol[ms2hsgenes$HGNC.symbol %in% common_genes])
length(mouse_common_genes)

seu_sc <- seu_sc[mouse_common_genes, ]
seuList <- lapply(seuList, function(seu) {
    seu[common_genes, ]
})

### 1.3.2 Normalize
seu_sc <- seu_sc %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000)
seuList <- lapply(seuList, function(seu_sp) {
    seu_sp %>%
        NormalizeData() %>%
        FindVariableFeatures(nfeatures = 2000)
})

### 1.3.3 align variable genes
mouse_vf <- VariableFeatures(seu_sc)
human_vf <- lapply(seuList, VariableFeatures)
h2m_vf <- lapply(human_vf, function(vf) {
    unique(ms2hsgenes$MGI.symbol[ms2hsgenes$HGNC.symbol %in% vf])
})

m_vf <- intersect(mouse_vf, Reduce(union, h2m_vf))
print(length(m_vf))

h_vf <- intersect(
    unique(ms2hsgenes$HGNC.symbol[ms2hsgenes$MGI.symbol %in% m_vf]),
    Reduce(union, human_vf)
)
print(length(h_vf))

VariableFeatures(seu_sc) <- m_vf
seuList <- lapply(seuList, function(seu_sp) {
    VariableFeatures(seu_sp) <- h_vf
    seu_sp
})

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
    marker <- marker.select(sig, overlap.max = 1)
    marker
})

n.marker <- sapply(markerList, function(marker) length(marker[[1]]))
message("The number of markers for each reference dataset are (", paste0(n.marker, collapse = ", "), ")")

save(markerList, file = "markerList.rda")

################################################################################
# 3. save spatial data
################################################################################
save(seuList, file = "seuList_processed.rda")

print(sapply(seuList, function(seu) apply(Embeddings(seu, "pos"), 2, max) - apply(Embeddings(seu, "pos"), 2, min)))

################################################################################
# 4. enrichment
################################################################################
library(msigdbr)

gene.used <- h_vf

pathwayList <- lapply(1:6, function(j) {
    hm_KEGG <- switch(j,
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP"),
        msigdbr(species = "Homo sapiens", category = "C4", subcategory = "CM"),
        msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
        msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB")
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
