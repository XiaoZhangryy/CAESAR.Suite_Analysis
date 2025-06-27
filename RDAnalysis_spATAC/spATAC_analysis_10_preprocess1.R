rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
setwd(work_wd)
source("functions.R")


figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)
enrich_wd <- paste0(work_wd, "Enrichment/")
dir.create(enrich_wd, showWarnings = FALSE)

setwd(work_wd)


################################################################################
# 1. QC scRNA-seq
################################################################################
MOCA.dir <- "/share/rawdata/zhangx/scRNA/MOCA/"
data_wd <- "/share/analysisdata/zhangx/spATAC/scRNAData/"

meta.data.RNA <- read.csv(file = paste0(MOCA.dir, "cell_annotate.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- read.csv(file = paste0(MOCA.dir, "gene_annotate.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- gene.ANN.RNA[, "gene_short_name", drop = FALSE]

cds <- readRDS(paste0(MOCA.dir, "gene_count_cleaned_sampled_100k.RDS"))

MOCA <- CreateSeuratObject(counts = cds, project = "MOCA")
meta.data.RNA <- meta.data.RNA[colnames(MOCA), ]
meta.data.RNA <- meta.data.RNA[, c("Main_cell_type", "development_stage")]

MOCA <- AddMetaData(object = MOCA, metadata = meta.data.RNA)
stages <- unique(meta.data.RNA$development_stage)

print(stages)
# [1] 12.5 13.5 11.5 10.5  9.5

print(dim(MOCA))
# [1]  26183 100000

stage <- 11.5

MOCA_sub <- subset(MOCA, development_stage == stage)
MOCA_sub.raw.data <- as.matrix(GetAssayData(MOCA_sub, slot = "counts"))
MOCA_sub.raw.data <- as.data.frame(MOCA_sub.raw.data)
MOCA_sub.raw.data <- merge(gene.ANN.RNA, MOCA_sub.raw.data, by = 0, all = TRUE)

# remove duplicated genes
tt <- table(MOCA_sub.raw.data$gene_short_name)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x) {
    rows <- which(MOCA_sub.raw.data$gene_short_name == x)
    return(rows[2:length(rows)])
}
row_del <- unlist(lapply(name_rep, row_del_fun))
MOCA_sub.raw.data <- MOCA_sub.raw.data[-row_del, ]

row.names(MOCA_sub.raw.data) <- MOCA_sub.raw.data$gene_short_name
MOCA_sub.raw.data <- MOCA_sub.raw.data[, -c(1:2), drop = FALSE]

print(dim(MOCA_sub.raw.data))
# [1] 26158 33802

seu_sc <- CreateSeuratObject(
    counts = MOCA_sub.raw.data,
    project = "MOCA_E11",
    meta.data = MOCA_sub@meta.data,
    min.cells = 1,
    min.features = 5
)

print(dim(seu_sc))
# [1] 25758 33802

rm(MOCA_sub.raw.data)

Idents(seu_sc) <- seu_sc$Main_cell_type
seu_sc <- subset(
    seu_sc,
    idents = c(
        "Lens", "Jaw and tooth progenitors", "Cardiac muscle lineages",
        "Melanocytes",
        "Myocytes", "Schwann cell precursor", # unrelated cell types
        "Cholinergic neurons", "Megakaryocytes" # ratio < 5e-3
    ), invert = TRUE
)
seu_sc$merged_clusters <- Idents(seu_sc)

print(dim(seu_sc))

save(seu_sc, file = "seu_sc_raw.rda")






