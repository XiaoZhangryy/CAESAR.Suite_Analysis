rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMxSimu/"
setwd(work_wd)
source("functions.R")

sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)

i <- commandArgs(TRUE) %>% as.integer()
ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
setwd(ana_wd)
load(paste0(ana_wd, "seu_all_raw.rda"))

seu <- seu %>%
    subset(downsample = 2000) %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 960)

n_cts <- table(Idents(seu))

################################################################################
# markers - generate once and keep
################################################################################
# deg_list <- lapply(1:4, function(i) {
#     load(paste0(work_wd, "processed_data_CosMx", i, "/degs.rda"))
#     degs
# })

# markerList <- find_marker(deg_list, 5, 1)

# save(markerList, file = "markerList.rda")

load("markerList.rda")

################################################################################
# run co-embeddings
################################################################################
## caesar
img_meta_wd <- paste0(ana_wd, "img_features/")
feature_img <- Reduce(rbind, lapply(fovs, function(fov) {
    meta_data <- read.csv(paste0(img_meta_wd, "Lung_meta", fov, ".csv"))
    feature_img <- fread(paste0(img_meta_wd, "feature_img_slice", fov, ".csv"))
    feature_img <- as.matrix(feature_img)
    rownames(feature_img) <- meta_data$X
    feature_img
}))
feature_img <- feature_img[colnames(seu), ]
save(feature_img, file = "feature_img.rda")

load("feature_img.rda")
feature_img <- feature_img[colnames(seu), ]
pos <- Embeddings(seu, "pos")
seu <- CAESAR.coembedding.image(
    seu, feature_img, pos,
    q = 10, reduction.name = "caesar", radius.upper = 200
)

## mca
seu <- RunMCA(seu, nmcs = 10, reduction.name = "mca")

## gedensity
ce <- compute.mca(object = seu, dims.use = 1:10)
el <- compute.nn.edges(coembed = ce)


################################################################################
# enrichment scoring
################################################################################
## cellid
pathway_cellid <- RunCellHGT(seu, pathways = markerList, dims = 1:10, minSize = 2, reduction = "mca")

## gsdenisty
pathway_gsdenisty <- run.rwr.list(el = el, gene_set_list = markerList, cells = colnames(seu))

## caesar
pathway_scores <- CAESAR.enrich.score(seu, pathwaylist = markerList)

pathway_gsea <- GSEA_enrich(seu, pathways = markerList)

pathway_gsva <- GSVA_enrich(seu, pathways = markerList)

pathway_aucell <- AUCell_enrich(seu, pathways = markerList)

pathway_vam <- VAM_enrich(seu, pathways = markerList)




res_auc <- Reduce(cbind, list(
    "CAESAR" = calculate_auc(seu$merged_clusters, pathway_scores, return.mean = FALSE),
    "Cell-ID" = calculate_auc(seu$merged_clusters, as.matrix(t(pathway_cellid)), return.mean = FALSE),
    "GSDensity" = calculate_auc(seu$merged_clusters, as.matrix(pathway_gsdenisty), return.mean = FALSE),
    "GSEA" = calculate_auc(seu$merged_clusters, pathway_gsea, return.mean = FALSE),
    "GSVA" = calculate_auc(seu$merged_clusters, pathway_gsva, return.mean = FALSE),
    "AUCell" = calculate_auc(seu$merged_clusters, pathway_aucell, return.mean = FALSE),
    "VAM" = calculate_auc(seu$merged_clusters, pathway_vam, return.mean = FALSE)
))
colnames(res_auc) <- c("CAESAR", "Cell-ID", "GSDensity")

save(res_auc, n_cts, file = paste0(ana_wd, "enrich_measures_all.rda"))


