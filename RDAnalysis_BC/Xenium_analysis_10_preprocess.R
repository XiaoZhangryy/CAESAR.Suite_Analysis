rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
dir.create(work_wd, showWarnings = FALSE)
setwd(work_wd)

source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

cols <- setNames(
    c(
        "#fdc086", "#386cb0", "#b30000", "#FBEA2E", "#731A73",
        "#FF8C00", "#F898CB", "#4DAF4A", "#a6cee3", "#737373"
    ),
    c(
        "B-cells", "CAFs", "Cancer Epithelial", "Endothelial", "Myeloid",
        "Normal Epithelial", "Plasmablasts", "PVL", "T-cells", "unassigned"
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
read_scRNA_BC <- function() {
    wd_here <- getwd()
    setwd("/share/rawdata/liuw/Projects/IntTemporalSpatial/RealData/xenium_hBreast2/Reference_scRNA/Wu_etal_2021_BRCA_scRNASeq/")

    expdat <- ReadMtx(
        "count_matrix_sparse.mtx",
        cells = "count_matrix_barcodes.tsv",
        features = "count_matrix_genes.tsv", feature.column = 1
    )
    metadat <- fread("metadata.csv")
    rowname_meta <- metadat$V1
    metadat <- as.data.frame(metadat)[, -1]
    rownames(metadat) <- rowname_meta

    print(dim(expdat))

    seu_sc <- CreateSeuratObject(
        counts = expdat, meta.data = metadat, min.features = 5, min.cells = 1
    )
    Idents(seu_sc) <- seu_sc$celltype_major
    seu_sc$merged_clusters <- Idents(seu_sc)

    print(dim(seu_sc))

    setwd(wd_here)

    return(seu_sc)
}

seu_sc <- read_scRNA_BC()

table(seu_sc$merged_clusters)

## -----------------------------------------------------------------------------
## 1.2 read Xenium data
## -----------------------------------------------------------------------------
read_sp_BC <- function() {
    load("/share/analysisdata/zhangx/Xenium_CoFAST/before_240606/AnalysisData/xenium_BC4_seuList.rds")

    print(sapply(seuList, dim))

    seuList <- lapply(seuList, function(seu) {
        seu <- CreateSeuratObject(
            counts = seu@assays$RNA@counts, meta.data = seu@meta.data,
            min.features = 5, min.cells = 1
        )

        cords <- c("x_centroid", "y_centroid")
        pos <- as.matrix(seu@meta.data[, cords])
        colnames(pos) <- paste0("pos", 1:2)

        seu@reductions[["pos"]] <- CreateDimReducObject(
            embeddings = pos,
            key = paste0("pos", "_"), assay = "RNA"
        )
        seu
    })

    seuList
}

seuList <- read_sp_BC()

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
VariableFeatures(seu_sc) <- common_genes
seuList <- lapply(seuList, function(seu_sp) {
    VariableFeatures(seu_sp) <- common_genes
    seu_sp
})


################################################################################
# 2. extract signature genes on reference
################################################################################
## 2.1 co-embedding
seu_sc <- NCFM(seu_sc, q = q_est)

save(seu_sc, file = "seu_sc.rda")

## 2.2 split seu_sc, rename, co-embedding, find sg, find marker
### 2.2.1 split seu_sc
seu_scList <- SplitObject(seu_sc, split.by = "orig.ident")

### 2.2.2 rename
sc_ct_ratios <- sapply(seu_scList, function(seu) {
    counts <- table(seu$merged_clusters)
    ratio <- setNames(rep(0, length(celltypes)), celltypes)
    ratio[names(counts)] <- counts / sum(counts)
    ratio
})
sc_ct_ratios <- sc_ct_ratios[names(sort(rowSums(sc_ct_ratios), decreasing = TRUE)), ]

reorder_matrix <- function(matrix, cell_types) {
    if (length(cell_types) != 2) {
        stop("Cell types vector must be of length 2.")
    }

    if (!all(cell_types %in% rownames(matrix))) {
        stop("Some cell types in the vector are not present in the matrix row names.")
    }

    order_columns <- function(matrix, cell_types) {
        primary_row_values <- matrix[cell_types[1], ]
        secondary_row_values <- matrix[cell_types[2], ]

        order <- order(-primary_row_values, secondary_row_values)

        return(order)
    }

    new_order <- order_columns(matrix, cell_types)
    reordered_matrix <- matrix[, new_order]

    return(reordered_matrix)
}

sc_ct_ratios <- reorder_matrix(sc_ct_ratios, c("Cancer Epithelial", "T-cells"))
sc_ct_ratios[1:2, ] <- sc_ct_ratios[2:1, ]
rownames(sc_ct_ratios)[1:2] <- rownames(sc_ct_ratios)[2:1]

sc_sample_id <- setNames(paste0("Sample ", 1:26), colnames(sc_ct_ratios))
names(seu_scList) <- sc_sample_id[names(seu_scList)]
seu_scList <- seu_scList[paste0("Sample ", 1:26)]

### 2.2.3 split co-embedding
seu_scList <- lapply(seu_scList, function(seu) {
    sample_id <- sc_sample_id[seu$orig.ident[1]] %>% as.character()
    seu$sample <- sample_id

    seu <- NCFM(seu, q = q_est)

    seu
})
save(seu_scList, file = "seu_scList.rda")

n_sc <- sapply(seu_scList, ncol)
print(n_sc)

print(median(n_sc))

### 2.2.4 find sg
sg_sc_List <- lapply(seu_scList, function(seu) {
    seu <- pdistance(seu, reduction = "ncfm")
    find.sig.genes(seu)
})
save(sg_sc_List, file = "sg_sc_List.rda")

### 2.2.5 find marker
markerList <- lapply(sg_sc_List, function(sig) {
    marker <- marker.select(sig, overlap.max = 1)
    if (length(marker[[1]]) <= 1) {
        marker <- marker.select(sig, overlap.max = 2)
    }
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

gene.used <- VariableFeatures(seu_sc)

pathwayList <- lapply(1:5, function(j) {
    hm_KEGG <- switch(j,
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
        msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP"),
        msigdbr(species = "Homo sapiens", category = "C4", subcategory = "CM"),
        msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
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




