

# ==============================================================================
# reference
# ArchRtoSignac: https://github.com/swaruplabUCI/ArchRtoSignac
# ArchR: https://www.archrproject.com/
# Spatial_ATAC-seq: https://github.com/dyxmvp/Spatial_ATAC-seq
# Signac: https://stuartlab.org/signac/
# MOCA: https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
# ==============================================================================

library(Seurat)
library(hdf5r)
library(Matrix)
# library(Signac)
library(data.table)
library(ArchR)

data.dir <- "/share/rawdata/zhangx/SpatialATAC/GSE171943/GSM5238385_ME11_50um/"
ana.dir <- "/share/analysisdata/zhangx/CAESAR/SpATACME/spATAC/"
dir.create(ana.dir, showWarnings = FALSE)
setwd(ana.dir)

threads <- 16
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- paste0(data.dir, "GSM5238385_ME11_50um.fragments.tsv.gz")
sampleNames <- 'ME11'

## Create ArchRProject
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sampleNames,
    minTSS = 0,
    minFrags = 0,
    maxFrags = 1e+07,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    offsetPlus = 0,
    offsetMinus = 0,
    TileMatParams = list(tileSize = 5000)
)

ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = sampleNames,
  copyArrows = TRUE
)

proj


## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
filter.matrix = TRUE
image <- Read10X_Image(
    image.dir = paste0(data.dir, "spatial"),
    filter.matrix = filter.matrix)
image.coor <- paste0(sampleNames, "#", row.names(image@coordinates), "-1")
meta.data.spatial <- meta.data[image.coor, ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue


## unsupervised clustering
proj_in_tissue <- addIterativeLSI(
    ArchRProj = proj_in_tissue,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_unsupervised",
  resolution = 0.5,
  force = TRUE
)

proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p_cluster <- plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters_unsupervised", embedding = "UMAP", size = 1.5)
ggsave(
    file = "umap_Clusters_unsupervised_IterativeLSI.png",
    plot = p_cluster, width = 9, height = 5.5, units = "in", dpi = 200)

proj_in_tissue <- addImputeWeights(proj_in_tissue)

## Identify the marker genes for each cluster 
# markersGS <- getMarkerFeatures(
#   ArchRProj = proj_in_tissue, 
#   useMatrix = "GeneScoreMatrix", 
#   groupBy = "Clusters_unsupervised",
#   testMethod = "wilcoxon"
# )

# markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

# markerGenes <- list()
# for (i in seq_len(length(markerList_pos))) {
#   markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
# }
# markerGenes <- unlist(markerGenes)

# save(markerGenes, file = "markerGenes_GeneScoreMatrix.rda")


## Call peaks
# proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters_unsupervised")

# pathToMacs2 <- findMacs2()

# proj_in_tissue <- addReproduciblePeakSet(
#   ArchRProj = proj_in_tissue,
#   groupBy = "Clusters",
#   pathToMacs2 = pathToMacs2,
#   force = TRUE
# )

# proj_in_tissue <- addPeakMatrix(proj_in_tissue)

# if ("Motif" %ni% names(proj_in_tissue@peakAnnotation)) {
#   proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
# }


slotNames(proj_in_tissue)

dim(proj_in_tissue@cellColData)

colnames(proj_in_tissue@cellColData)

names(proj_in_tissue@reducedDims)

getAvailableMatrices(proj_in_tissue)


## load scRNA data - MOCA data
load("/share/analysisdata/zhangx/CAESAR/SpATACME/seu_sc_raw.rda")

# ------------------------------------------------------------------------------
## Integration with ArchR oject
# ------------------------------------------------------------------------------
proj_in_tissue <- addGeneIntegrationMatrix(
    ArchRProj = proj_in_tissue, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seu_sc,
    addToArrow = TRUE,
    groupRNA = "Main_cell_type",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    force = TRUE
)

# proj_in_tissue <- addGeneIntegrationMatrix(
#     ArchRProj = proj_in_tissue, 
#     useMatrix = "GeneScoreMatrix",
#     matrixName = "GeneIntegrationMatrix2",
#     reducedDims = "IterativeLSI",
#     seRNA = seu_sc,
#     addToArrow = TRUE,
#     groupRNA = "merged_clusters",
#     nameCell = "predictedCell2",
#     nameGroup = "predictedGroup2",
#     nameScore = "predictedScore2",
#     force = TRUE
# )




meta <- data.frame(proj_in_tissue@cellColData)
colnames(meta)
new_row_names <- row.names(meta)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
rownames(meta) <- new_row_names

positions <- read.csv(file = paste0(data.dir, 'spatial/tissue_positions_list.csv'), header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(positions) <- c("in_tissue", "array_row", "array_column", "pxl_col_in_fullres", "pxl_row_in_fullres")
positions <- positions[new_row_names, ]

meta <- cbind(meta, positions)

## build ATAC seurat object
useMatrix <- "GeneScoreMatrix"
reducedDims <- "IterativeLSI"
dimsToUse <- 1:30
scaleDims <- NULL
corCutOff <- 0.75

subProj <- proj_in_tissue
subProj@imputeWeights <- SimpleList()
geneDF <- ArchR:::.getFeatureDF(getArrowFiles(subProj), useMatrix)
# > dim(geneDF)
# [1] 24333     6
geneDF <- geneDF[geneDF$name %in% rownames(seu_sc), , drop = FALSE]
# > dim(geneDF)
# [1] 18328     6

mat <- ArchR:::.getPartialMatrix(
    getArrowFiles(subProj),
    featureDF = geneDF,
    threads = 1,
    cellNames = subProj$cellNames,
    useMatrix = useMatrix,
    verbose = FALSE
)
rownames(mat) <- geneDF[, "name"]

lsi_embed <- getReducedDims(ArchRProj = subProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
rownames(lsi_embed) <- new_row_names

# mat <- log(mat + 1) #use natural log
mat0 <- mat
colnames(mat) <- new_row_names
seuATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
seuATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
DefaultAssay(seuATAC) <- "GeneScore"
seuATAC@assays$GeneScore@data <- as.matrix(log(mat + 1))

seuATAC[["lsi"]] <- CreateDimReducObject(
    embeddings = lsi_embed,
    key = "LSI_",
    assay = "GeneScore"
)
seuATAC <- AddMetaData(seuATAC, meta)

print(seuATAC)

save(seuATAC, file = "seuATAC_final.rda")

saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-ProjSpATAC", load = FALSE)
