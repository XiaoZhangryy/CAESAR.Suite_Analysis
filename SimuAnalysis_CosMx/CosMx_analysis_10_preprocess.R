rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMxSimu/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)
source("functions.R")

################################################################################
# 1. preprocess
################################################################################
i <- commandArgs(TRUE) %>% as.integer()
sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)
sample <- sample_names[i]

ana_folder <- paste0("processed_data_CosMx", i)
ana_wd <- paste0(work_wd, ana_folder, "/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

## sample i
load(paste0(
    "/share/rawdata/zhangx/CosMx/HumanNSCLC/processed_data_Seurat/Seurat_objs/",
    sample, "_Seurat_obj.rda"
))
seu <- subset(seu, subset = nFeature_RNA >= 5 & nFeature_RNA < 500)
Idents(seu) <- seu$cell_type

newcluster.ids <- c(
    "B-cell", "endothelial", "epithelial", "fibroblast", "macrophage", "mast",
    "mDC", "monocyte", "neutrophil", "NK", "pDC", "plasmablast", "T-cell",
    "T-cell", "T-cell", "T-cell", "T-cell", "tumor",
    "tumor", "tumor", "tumor", "tumor"
)
names(newcluster.ids) <- sort(levels(seu))
seu <- RenameIdents(seu, newcluster.ids)
seu$merged_clusters <- Idents(seu)

cords <- c("CenterX_global_px", "CenterY_global_px")
pos <- as.matrix(seu@meta.data[, cords])
colnames(pos) <- paste0("pos", 1:2)

seu@reductions[["pos"]] <- CreateDimReducObject(
    embeddings = pos,
    key = paste0("pos", "_"), assay = "RNA"
)

save(seu, file = "seu_all_raw.rda")

################################################################################
# 2. degs
################################################################################
seu <- seu %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 960)

degs <- FindAllMarkers(seu)
save(degs, file = "degs.rda")


################################################################################
# 3. degs for bin 5 data
################################################################################
bin5pos <- segment_square(pos, sq_nspots = sqrt(ncol(seu) / 5) + 1, by_order = TRUE, verbose = TRUE)
seubin5 <- get_merged_seu(seu, bin5pos$spotID_list, bin5pos$pos_new)
Idents(seubin5) <- seubin5$merged_clusters

seubin5 <- seubin5 %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 960)

degsbin5 <- FindAllMarkers(seubin5)
save(degsbin5, file = "degsbin5.rda")

################################################################################
# 4. Split by Field of View (FoV), co-embedding with CAESAR and MCA
################################################################################
celltypes <- sort(levels(seu))
cols <- setNames(
    PRECAST::chooseColors(
        palettes_name = "Classic 20", n_colors = length(celltypes)
    ),
    celltypes
)
save(cols, file = "cols.rda")
save(celltypes, file = "celltypes.rda")
q_est <- 50
save(q_est, file = "q_est.rda")


print(table(Idents(seu)))

## # Split data by FoV, Co-embedding with caesar and MCA
seuList <- SplitObject(seu, split.by = "fov")

fovs <- names(seuList)
save(fovs, file = "fovs.rda")

img_meta_wd <- paste0(ana_wd, "img_features/")
for (fov in fovs) {
    fov_wd <- paste0(ana_wd, "fov", fov, "/")
    dir.create(fov_wd, showWarnings = FALSE)
    setwd(fov_wd)

    seu_i <- seuList[[fov]]
    pos <- Embeddings(seu_i, "pos")
    meta_data <- read.csv(paste0(img_meta_wd, "Lung_meta", fov, ".csv"))
    set.seed(2024)
    feature_img <- fread(paste0(img_meta_wd, "feature_img_slice", fov, ".csv"))
    feature_img <- as.matrix(feature_img)
    rownames(feature_img) <- meta_data$X

    # Run CAESAR co-embedding
    seu_i <- tryCatch(
        {
            CAESAR.coembedding.image(
                seu_i, feature_img, pos,
                q = q_est, reduction.name = "caesar", radius.upper = 100
            )
        },
        error = function(cond) {
            CAESAR.coembedding.image(
                seu_i, feature_img, pos,
                q = q_est, reduction.name = "caesar", radius.upper = 150
            )
        }
    )
    seu_i <- pdistance(seu_i, reduction = "caesar")

    sig_list <- find.sig.genes(seu_i)
    sig_list <- list(sig_list)
    save(sig_list, file = paste0("sig_list_scRNA_fov", fov, ".rda"))
    celltypes <- sort(levels(seu_i))
    save(celltypes, file = paste0("celltypes_fov", fov, ".rda"))

    # RunMCA
    seu_i <- RunMCA(seu_i, nmcs = q_est, reduction.name = "mca")

    # Save original FoV-level Seurat object
    seu <- seu_i
    seu@assays[["distce"]] <- NULL
    seu_name <- paste0("seu_fov", fov, ".rda")
    save(seu, file = seu_name)

    message("sample fov ", fov, " done.")




    # Repeat CAESAR pipeline for bin5 version
    bin5pos <- segment_square(pos, sq_nspots = sqrt(ncol(seu) / 5) + 1, by_order = TRUE, verbose = TRUE)
    merged_res <- get_merged_seu_fov(seu, bin5pos$spotID_list, bin5pos$pos_new, feature_img)
    seubin5 <- merged_res$seu
    feature_imgbin5 <- merged_res$feature_img
    posbin5 <- as.matrix(seubin5@meta.data[, c("row", "col")])

    colnames(posbin5) <- paste0("pos", 1:2)
    seubin5@reductions[["pos"]] <- CreateDimReducObject(
        embeddings = posbin5,
        key = paste0("pos", "_"), assay = "RNA"
    )
    Idents(seubin5) <- seubin5@meta.data$merged_clusters

    seubin5 <- seubin5 %>%
        NormalizeData() %>%
        FindVariableFeatures(nfeatures = 960)

    seu_i <- CAESAR.coembedding.image(
        seubin5, feature_imgbin5, posbin5,
        q = q_est, reduction.name = "caesar", radius.upper = 400
    )
    seu_i <- pdistance(seu_i, reduction = "caesar")

    sig_list <- find.sig.genes(seu_i)
    sig_list <- list(sig_list)
    save(sig_list, file = paste0("sig_list_scRNA_fov", fov, "_bin5.rda"))
    celltypes <- sort(levels(seu_i))
    save(celltypes, file = paste0("celltypes_fov", fov, "_bin5.rda"))


    seu_i <- RunMCA(seu_i, nmcs = q_est, reduction.name = "mca")

    seu <- seu_i
    seu@assays[["distce"]] <- NULL
    seu_name <- paste0("seu_fov", fov, "_bin5.rda")
    save(seu, file = seu_name)

    message("sample fov ", fov, " bin5 done.")
}


 
