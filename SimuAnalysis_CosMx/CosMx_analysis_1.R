rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMx/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

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


################################################################################
# 2. co-embedding
################################################################################
seu <- seu %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 960)

degs <- FindAllMarkers(seu)
save(degs, file = "degs.rda")

load("feature_img.rda")
seu <- CAESAR.coembedding.image(
    seu, feature_img, pos,
    q = 50, reduction.name = "caesar", radius.upper = 100
)
save(seu, file = "seu_all.rda")

celltypes <- sort(levels(seu))
cols <- setNames(
    PRECAST::chooseColors(
        palettes_name = "Classic 20", n_colors = length(celltypes)
    ),
    celltypes
)
save(cols, file = "cols.rda")
save(celltypes, file = "celltypes.rda")
save(q_est, file = "q_est.rda")

q_est <- 50

print(table(Idents(seu)))

### 1.4.2 Co-embedding with caesar
seuList <- SplitObject(seu, split.by = "fov")

fovs <- names(seuList)
save(fovs, file = "fovs.rda")

load("feature_img_list.rda")
for (fov in fovs) {
    seu_i <- seuList[[fov]]
    pos <- Embeddings(seu_i, "pos")
    feature_img <- feature_img_list[[fov]]
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

    # cellid
    seu_i <- RunMCA(seu_i, nmcs = q_est, reduction.name = "mca")

    seu <- seu_i
    seu@assays[["distce"]] <- NULL
    seu_name <- paste0(ana_wd, "seu_fov", fov, ".rda")
    save(seu, file = seu_name)

    message("sample fov ", fov, " done.")
}





################################################################################
# 3. bin 5 data
################################################################################
segment_square <- function(pos, sq_nspots = 70, by_order = T, verbose = T) {
    tmp <- pos[, 1]
    if (by_order) {
        x_cut <- sort(tmp)[seq(1, length(tmp), length.out = sq_nspots + 1)]
    } else {
        x_cut <- seq(min(tmp), max(tmp), length.out = sq_nspots + 1)
    }

    tmp <- pos[, 2]
    if (by_order) {
        y_cut <- sort(tmp)[seq(1, length(tmp), length.out = sq_nspots + 1)]
    } else {
        y_cut <- seq(min(tmp), max(tmp), length.out = sq_nspots + 1)
    }

    i <- 1
    pos_new <- matrix(NA, sq_nspots^2, 2)
    areaList <- list()
    for (i1 in 1:sq_nspots) {
        if (verbose) {
            message("i1 = ", i1)
        }
        for (i2 in 1:sq_nspots) {
            if (i1 < sq_nspots && i2 < sq_nspots) {
                tmp <- which(x_cut[i1] <= pos[, 1] & pos[, 1] < x_cut[i1 + 1] & y_cut[i2] <= pos[, 2] & pos[, 2] < y_cut[i2 + 1])
            } else if (i1 < sq_nspots && i2 == sq_nspots) {
                tmp <- which(x_cut[i1] <= pos[, 1] & pos[, 1] < x_cut[i1 + 1] & y_cut[i2] <= pos[, 2] & pos[, 2] <= y_cut[i2 + 1])
            } else {
                tmp <- which(x_cut[i1] <= pos[, 1] & pos[, 1] <= x_cut[i1 + 1] & y_cut[i2] <= pos[, 2] & pos[, 2] < y_cut[i2 + 1])
            }

            areaList[[i]] <- tmp
            pos_new[i, ] <- c((x_cut[i1] + x_cut[i1 + 1]) / 2, (y_cut[i2] + y_cut[i2 + 1]) / 2)
            i <- i + 1
        }
    }
    idx <- which(sapply(areaList, function(x) length(x) > 0))
    return(list(spotID_list = areaList[idx], pos_new = pos_new[idx, ]))
}

get_merged_seu <- function(seu, areaList, pos_new) {
    require(Seurat)
    n_area <- length(areaList)
    count_new <- matrix(NA, nrow(seu), n_area)
    colnames(count_new) <- paste0("merge_spot", 1:n_area)
    row.names(count_new) <- row.names(seu)
    merged_clusters <- rep(NA, n_area)
    DefaultAssay(seu) <- "RNA"
    for (i in 1:n_area) { #
        message("i = ", i, "/", n_area)
        if (length(areaList[[i]]) > 1) {
            count_new[, i] <- rowSums(seu[["RNA"]]@counts[, areaList[[i]]])
        } else {
            count_new[, i] <- seu[["RNA"]]@counts[, areaList[[i]]]
        }
        count_ct <- table(seu$merged_clusters[areaList[[i]]])
        domains <- names(which(count_ct == max(count_ct)))
        merged_clusters[i] <- domains[sample(seq_along(domains), 1)]
    }
    rm(seu)
    meta_data <- data.frame(
        row = pos_new[, 1], col = pos_new[, 2],
        merged_clusters = merged_clusters
    )
    row.names(meta_data) <- colnames(count_new)
    seu_new <- CreateSeuratObject(counts = as.sparse(count_new), meta.data = meta_data)
    seu_new
}

bin5pos <- segment_square(pos, sq_nspots = sqrt(ncol(seu) / 5) + 1, by_order = TRUE, verbose = TRUE)
seubin5 <- get_merged_seu(seu, bin5pos$spotID_list, bin5pos$pos_new)
Idents(seubin5) <- seubin5$merged_clusters

seubin5 <- seubin5 %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 960)

degsbin5 <- FindAllMarkers(seubin5)
save(degsbin5, file = "degsbin5.rda")





img_meta_wd <- paste0("/share/analysisdata/liuw/coembed_image/RealData/LungCancerCosMix/processedMeta/sample", i, "/")
for (fov in fovs) {
    load(paste0(ana_wd, "seu_fov", fov, ".rda"))
    pos <- Embeddings(seu, "pos")
    meta_data <- read.csv(paste0(img_meta_wd, "Lung_meta", fov, ".csv"))
    message("bardcode fro fov ", fov, " aligned: ", all(meta_data$barcode == colnames(seu)))
    set.seed(2024)
    feature_img <- fread(paste0(img_meta_wd, "feature_img_slice", fov, ".csv"))
    feature_img <- as.matrix(feature_img)
    rownames(feature_img) <- meta_data$barcode

    bin5pos <- segment_square(pos, sq_nspots = sqrt(ncol(seu) / 5) + 1, by_order = TRUE, verbose = TRUE)
    merged_res <- get_merged_seu(seu, bin5pos$spotID_list, bin5pos$pos_new, feature_img)
    seubin5 <- merged_res$seu
    feature_imgbin5 <- merged_res$feature_img
    posbin5 <- as.matrix(seubin5@meta.data[, c("row", "col")])

    colnames(posbin5) <- paste0("pos", 1:2)
    seubin5@reductions[["pos"]] <- CreateDimReducObject(
        embeddings = posbin5,
        key = paste0("pos", "_"), assay = "RNA"
    )
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
    seu_name <- paste0(ana_wd, "seu_fov", fov, "_bin5.rda")
    save(seu, file = seu_name)

    message("sample fov ", fov, " done.")
}
