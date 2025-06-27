rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/STMOB_coarse/"
setwd(work_wd)
source("functions.R")

dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)

cols_manual <- setNames(
    c(
        "#4374A5", "#FCDDDE", "#2AB67F", "#F08A21", "#737373"
    ),
    c(
        "GCL", "MCL", "ONL", "GL", "Unknown"
    )
)

cols <- setNames(
    c(
        "#4374A5", "#FCDDDE", "#2AB673", "#F08A21", "#E04D50",
        "#737373"
    ),
    c(
        "GC", "M/TC", "OSNs", "PGC", "EPL-IN", "unassigned"
    )
)
celltypes <- setdiff(names(cols), "unassigned")
q_est <- 50
save(cols, cols_manual, celltypes, q_est, file = "cols_cts_q.rda")

################################################################################
# 1. preprocess
################################################################################
## -----------------------------------------------------------------------------
## 1.1 read scRNA data
## -----------------------------------------------------------------------------
read_scRNA_STMOB <- function() {
    load("/share/analysisdata/zhangx/MOB_CoFAST/AnalysisData/MOB.dge.sceset.RData")

    counts <- assay(sce, "counts")
    meta.data <- as.data.frame(colData(sce))

    print(dim(counts))

    seu_sc <- CreateSeuratObject(
        counts = counts, meta.data = meta.data,
        min.features = 5, min.cells = 1
    )
    Idents(seu_sc) <- seu_sc@meta.data$cellType
    seu_sc$merged_clusters <- Idents(seu_sc)

    print(dim(seu_sc))

    return(seu_sc)
}

seu_sc <- read_scRNA_STMOB()

table(seu_sc$merged_clusters)

## -----------------------------------------------------------------------------
## 1.2 read Xenium data
## -----------------------------------------------------------------------------
read_sp_STMOB <- function() {
    load("/share/analysisdata/zhangx/MOB_CoFAST/AnalysisData/Rep12_MOB_count_matrix-1.RData")
    load("/share/analysisdata/zhangx/MOB_CoFAST/AnalysisData/Rep12_MOB_manual_annotation.RData")

    coordinate <- read.csv("/share/analysisdata/zhangx/MOB_CoFAST/AnalysisData/coordinate.csv", header = FALSE)
    meta.data <- cbind(coordinate, Rep12_MOB_manual_annotation)
    colnames(meta.data) <- c("x", "y", "manual_annotation")
    rownames(meta.data) <- colnames(MOB_raw)

    print(dim(MOB_raw))

    seu <- CreateSeuratObject(
        counts = MOB_raw, meta.data = meta.data,
        min.features = 1, min.cells = 1
    )
    Idents(seu) <- seu@meta.data$manual_annotation

    cords <- c("x", "y")
    pos <- as.matrix(seu@meta.data[, cords])
    colnames(pos) <- paste0("pos", 1:2)

    seu@reductions[["pos"]] <- CreateDimReducObject(
        embeddings = pos,
        key = paste0("pos", "_"), assay = "RNA"
    )

    print(dim(seu))

    seu
}

seu <- read_sp_STMOB()

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
    FindVariableFeatures(nfeatures = 2000)
seu <- seu %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000)

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
expr.prop.cutoff <- 0.1
markerList <- lapply(sg_sc_List, function(sig) {
    marker <- marker.select(sig, expr.prop.cutoff = expr.prop.cutoff, ntop.max = 200, overlap.max = 1)
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
