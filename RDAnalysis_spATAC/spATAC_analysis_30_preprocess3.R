rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
enrich_wd <- paste0(work_wd, "Enrichment/")

ana_wd <- paste0(work_wd, "spATACME/")
dir.create(ana_wd, showWarnings = FALSE)

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
## 1.1 align scRNA and Xenium data
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
expr.prop.cutoff <- 0.1
markerList <- lapply(sg_sc_List, function(sig) {
    marker <- marker.select(sig, expr.prop.cutoff = expr.prop.cutoff, ntop.max = 200, overlap.max = 4)
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




save(seu, file = paste0(ana_wd, "seu.rda"))

seu_sc@assays[["distce"]] <- NULL
save(seu_sc, file = paste0(ana_wd, "seu_sc.rda"))
