rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/HmHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

cols <- setNames(
    c(
        "#fdc086", "#386cb0", "#F08A21", "#b30000", "#70B5B0",
        "#a6cee3", "#ffff99", "#737373"
    ),
    c(
        "B cell", "CAF", "HPC-like", "Malignant cell",
        "T cell", "TAM", "TEC", "unassigned"
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
read_scRNA_hmHCC <- function() {
    counts <- Read10X(data.dir = "/share/analysisdata/zhangx/HCC_CoFAST/rawData")
    meta_sc <- read.csv("/share/analysisdata/zhangx/HCC_CoFAST/rawData/GSE125449_Set1_samples.txt", skip = 1, sep = "\t", header = FALSE)
    rownames(meta_sc) <- meta_sc$V2
    meta_sc$V2 <- NULL
    colnames(meta_sc) <- c("fov", "celltype")

    print(dim(counts))

    seu_sc <- CreateSeuratObject(
        counts = counts, meta.data = meta_sc,
        min.features = 5, min.cells = 1
    )
    Idents(seu_sc) <- seu_sc$celltype
    seu_sc$merged_clusters <- Idents(seu_sc)

    seu_sc <- seu_sc[, Idents(seu_sc) != "unclassified"]

    print(dim(seu_sc))

    return(seu_sc)
}

seu_sc <- read_scRNA_hmHCC()

table(seu_sc$merged_clusters)

## -----------------------------------------------------------------------------
## 1.2 read Xenium data
## -----------------------------------------------------------------------------
read_sp_hmHCC <- function() {
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

seuList <- read_sp_hmHCC()

## -----------------------------------------------------------------------------
## 1.3 align scRNA and Xenium data
## -----------------------------------------------------------------------------
common_genes <- Reduce(intersect, c(
    list(rownames(seu_sc)),
    lapply(seuList, rownames)
))
print(length(common_genes))

### 1.3.1 align genes
seu_sc <- seu_sc[common_genes, ]
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
common_vgs <- intersect(
    VariableFeatures(seu_sc),
    Reduce(union, lapply(seuList, VariableFeatures))
)
print(length(common_vgs))

VariableFeatures(seu_sc) <- common_vgs
seuList <- lapply(seuList, function(seu_sp) {
    VariableFeatures(seu_sp) <- common_vgs
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

gene.used <- common_vgs

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
