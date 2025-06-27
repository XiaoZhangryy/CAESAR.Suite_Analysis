library(Seurat)
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)
library(abind)

assign4CosMx <- function(number) {
    ns <- c(29, 28, 20)
    i <- 2
    while (length(ns) > 0) {
        if (number > ns[1]) {
            number <- number - ns[1]
            ns <- ns[-1]
            i <- i + 1
        } else {
            load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
            j <- as.integer(fovs[number])
            return(c(i, j))
        }
    }
}

Seurat.transfer <- function(seu, ref_seurat, true_y){
    set.seed(1)
    prediction <- tryCatch(
        {
            transfer_anchor <- FindTransferAnchors(reference = ref_seurat, query = seu)
            TransferData(anchorset = transfer_anchor, refdata = ref_seurat@meta.data[, true_y])
        }, error = function(msg){
            transfer_anchor <- FindTransferAnchors(reference = ref_seurat, query = seu, k.filter = NA)
            TransferData(anchorset = transfer_anchor, refdata = ref_seurat@meta.data[, true_y])
    })

    pred <- prediction$predicted.id
    significant <- prediction$prediction.score.max
    pred_un <- ifelse(significant > 0.5, pred, "unassigned")
  
    return(
        data.frame(
            pred = pred,
            pred_unassign = pred_un,
            significant = significant
        )
    )
}

scmap.transfer <- function(seu, ref_seurat, true_y) {
    library(scmap)
    set.seed(1)
    ref                                   <- as.SingleCellExperiment(ref_seurat)
    rowData(ref)$feature_symbol           <- rownames(ref)
    counts(ref)                           <- as.matrix(counts(ref))
    logcounts(ref)                        <- as.matrix(logcounts(ref))
    ref                                   <- selectFeatures(ref)
    ref                                   <- indexCluster(ref, cluster_col = true_y)
  index_cluster <- list((ref@metadata)$scmap_cluster_index)
    sce <- as.SingleCellExperiment(seu)
    rowData(sce)$feature_symbol           <- rownames(sce)
    counts(sce)                           <- as.matrix(counts(sce))
    logcounts(sce)                        <- as.matrix(logcounts(sce))
    scmap_cluster <- scmapCluster(sce, index_cluster, threshold = 0)

    pred <- scmap_cluster$combined_labs
    significant <- scmap_cluster$scmap_cluster_siml
    significant[is.na(significant)] <- -1
    pred_un <- ifelse(significant > 0.7, pred, "unassigned")

    return(
        data.frame(
            pred = pred,
            pred_unassign = pred_un,
            significant = significant
        )
    )
}

SingleR.transfer <- function(seu, ref_seurat, true_y) {
    library(SingleR)
    set.seed(1)
    sce <- as.SingleCellExperiment(ref_seurat)
    sce2 <- as.SingleCellExperiment(seu)
    genes <- intersect(rownames(sce), rownames(sce2))
    sce <- sce[genes, ]
    sce2 <- sce2[genes, ]
    SingleR_train <-
        trainSingleR(ref = sce, labels = colData(sce)[, true_y])
    pred_SingleR <-
        SingleR::classifySingleR(
            sce2,
            trained = SingleR_train
        )

    pred <- pred_SingleR$labels
    significant <- pred_SingleR$delta.next
    pred_un <- pred_SingleR$pruned.labels
    pred_un[is.na(pred_un)] <- "unassigned"

    return(
        data.frame(
            pred = pred,
            pred_unassign = pred_un,
            significant = significant
        )
    )
}

scPred.transfer <- function(x, ref_seurat, true_y) {
    library(scPred)
    set.seed(1)

    scp <- tryCatch(
        {
            ref_seurat %>% 
            ScaleData() %>%
            RunPCA() %>% 
            getFeatureSpace(true_y) %>%
            trainModel(seed = 1)
        }, error = function(msg){
            n_ct <- table(ref_seurat@meta.data[, true_y])
            cts <- names(n_ct)[n_ct > 5]
            ref_seurat <- ref_seurat[, ref_seurat@meta.data[, true_y] %in% cts]

            ref_seurat %>% 
            ScaleData() %>%
            RunPCA() %>% 
            getFeatureSpace(true_y) %>%
            trainModel(seed = 1)
    })
    x <- scPredict(x, scp)

    pred <- x$scpred_no_rejection
    significant <- x$scpred_max
    pred_un <- pred
    pred_un[significant < 0.55] <- "unassigned"

    return(
        data.frame(
            pred = pred,
            pred_unassign = pred_un,
            significant = significant
        )
    )
}

CelliD.Annotation <- function(seu, markerList) {
    q_est <- ncol(Embeddings(seu, "mca"))
    HGTres <- RunCellHGT(seu, pathways = markerList, dims = 1:q_est, minSize = 2, reduction = "mca")
    pred <- rownames(HGTres)[apply(HGTres, 2, which.max)]
    significant = apply(HGTres, 2, max)
    pred_un <- ifelse(significant >= 2, pred, "unassigned")

    return(
        data.frame(
            pred = pred,
            pred_unassign = pred_un,
            significant = significant
        )
    )
}

calculate_acc <- function(pred_y, true_y, cts_in_ref = NULL) {
    pred_y <- as.character(pred_y)
    true_y <- as.character(true_y)
    if (!is.null(cts_in_ref)) {
        id <- which(true_y %in% cts_in_ref)
    } else {
        id <- seq_along(pred_y)
    }
    return(mean(pred_y[id] == true_y[id]))
}

calculate_auc <- function(celltype, pathway.scores, return.mean = TRUE) {
    if (!is.factor(celltype) && !is.character(celltype)) {
        stop("'celltype' must be a factor or character vector representing cell type labels.")
    }

    if (!is.matrix(pathway.scores)) {
        stop("'pathway.scores' must be a matrix with cells as rows and pathways as columns.")
    }

    if (nrow(pathway.scores) != length(celltype)) {
        stop("The number of rows in 'pathway.scores' must match the length of 'celltype'.")
    }

    y <- as.character(celltype)
    cts <- colnames(pathway.scores)
    n <- length(y)

    res <- sapply(cts, function(ct) {
        n_ct <- sum(y == ct)
        if (n_ct == 0) {
            return(0)
        }

        xx <- pathway.scores[, ct]
        yy <- y == ct
        id <- order(xx, yy, decreasing = TRUE)

        x0 <- seq_len(n) / n
        y0 <- cumsum(yy[id]) / n_ct

        auc1 <- DescTools::AUC(x0, y0)

        id <- order(xx, -yy, decreasing = TRUE)
        y1 <- cumsum(yy[id]) / n_ct

        auc2 <- DescTools::AUC(x0, y1)

        return((auc1 + auc2) / 2)
    })

    if (return.mean) {
        weight <- table(y)[cts]
        id <- which((!is.na(weight)) & (weight != 0))
        weight <- weight[id] / sum(weight[id])
        return(as.vector(res[id] %*% weight))
    } else {
        return(res)
    }
}

accdiff <- function(pred_y, pred_y_unasg, true_y, cts_in_ref = NULL) {
    pred_y <- as.character(pred_y)
    pred_y_unasg <- as.character(pred_y_unasg)
    true_y <- as.character(true_y)

    if (!is.null(cts_in_ref)) {
        id <- which(true_y %in% cts_in_ref)
    } else {
        id <- seq_along(pred_y)
    }
    pred_y <- pred_y[id]
    pred_y_unasg <- pred_y_unasg[id]
    true_y <- true_y[id]

    id1 <- which(pred_y_unasg == "unassigned")

    acc1 <- mean(pred_y[id1] == true_y[id1])
    acc2 <- mean(pred_y == true_y)

    nn <- length(id1)

    res <- c(
        "n" = nn,
        "ACCunasg" = acc1,
        "ACC" = acc2,
        "ACCdiff" = acc2 - acc1
    )

    return(res)
}

asw <- function(embeddings, celltype) {
    n <- nrow(embeddings)
    if (n > 50000) {
        set.seed(1)
        id <- sample(1:n, 50000)
        embeddings <- embeddings[id, ]
        celltype <- celltype[id]
        rm(id)
    }
    metadata <- data.frame(celltype = celltype)
    sh_scores_pro <- scPOP::silhouette_width(embeddings, meta.data = metadata, c("celltype"))

    return(c(asw = sh_scores_pro[1]))
}

mymeasures <- function(pred_y, true_y, cts_in_ref) {
    pred_y <- as.character(pred_y)
    true_y <- as.character(true_y)

    POU <- mean(pred_y == "unassigned")
    ACC <- calculate_acc(pred_y, true_y, cts_in_ref)

    return(setNames(c(POU, ACC), c("POU", "ACC")))
}

calculate_f1_scores <- function(pred_y, true_y, cts_in_ref) {
    pred_y <- as.character(pred_y)
    true_y <- as.character(true_y)
    cm <- table(true_y, pred_y)

    cts_in_ref <- intersect(cts_in_ref, rownames(cm))
    cm <- cm[cts_in_ref, ]

    result <- sapply(cts_in_ref, function(x) {
        if (!(x %in% colnames(cm)) || ncol(cm) < 2) {
            return(c(Precision = 0, Recall = 0, F1 = 0))
        }
        false_xr <- setdiff(cts_in_ref, x)
        false_xc <- setdiff(colnames(cm), x)

        TP <- cm[x, x]
        FP <- sum(cm[false_xr, x])
        FN <- sum(cm[x, false_xc])

        Precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
        Recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
        F1 <- ifelse((Precision + Recall) == 0, 0, 2 * Precision * Recall / (Precision + Recall))

        return(c(Precision = Precision, Recall = Recall, F1 = F1))
    })

    n <- table(true_y)[cts_in_ref]

    return(rbind(
        result, n = n
    ))
}

pred_control_pou <- function(anno.df, pred_un, pred, pred_conf, pou_th = pou_th) {
    pred_y_un <- anno.df[, pred_un]
    pred_y_pou <- mean(pred_y_un == "unassigned")
    if (pred_y_pou < pou_th) {
        return(pred_y_un)
    }
    pred_y <- anno.df[, pred]
    pred_conf <- anno.df[, pred_conf]
    pred_y[pred_conf < quantile(pred_conf, pou_th)] <- "unassigned"
    return(pred_y)
}

find_marker <- function(deg_list, k = 5, log2fc = 1) {
    markerList_list <- lapply(deg_list, function(degs) {
        degs_top3 <- degs %>%
            filter(p_val_adj < 0.05, avg_log2FC > log2fc) %>%
            group_by(cluster)

        markerList <- split(degs_top3$gene, degs_top3$cluster)
        markerList
    })

    cell_types <- names(markerList_list[[1]])
    common_genes_list <- vector("list", length = length(cell_types))
    names(common_genes_list) <- cell_types
    for (ct in cell_types) {
        common_genes_list[[ct]] <- Reduce(c, lapply(markerList_list, function(xx) {
            xx[[ct]]
        }))

        freq_table <- table(common_genes_list[[ct]])
        freq_df <- as.data.frame(freq_table)
        colnames(freq_df) <- c("element", "frequency")
        freq_df <- freq_df[sample(nrow(freq_df)), ]
        freq_df <- freq_df[order(-freq_df$frequency), ]

        k0 <- freq_df$frequency[k]
        top_k_fix <- freq_df$element[freq_df$frequency > k0]
        top_k_candidate <- freq_df$element[freq_df$frequency == k0]

        if (length(top_k_fix) + length(top_k_candidate) == k) {
            common_genes_list[[ct]] <- c(top_k_fix, top_k_candidate)
        } else {
            top_k_candidate <- sample(top_k_candidate, k - length(top_k_fix))
            common_genes_list[[ct]] <- c(top_k_fix, top_k_candidate)
        }
    }

    common_genes_list <- common_genes_list[sapply(common_genes_list, length) > 0]

    return(common_genes_list)
}

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
    for (i in 1:n_area) { 
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

get_merged_seu_fov <- function(seu, areaList, pos_new, feature_img) {
    require(Seurat)
    n_area <- length(areaList)
    count_new <- matrix(NA, nrow(seu), n_area)
    colnames(count_new) <- paste0("merge_spot", 1:n_area)
    row.names(count_new) <- row.names(seu)
    merged_clusters <- rep(NA, n_area)
    feature_img_new <- matrix(NA, n_area, ncol(feature_img))
    DefaultAssay(seu) <- "RNA"
    for (i in 1:n_area) { 
        message("i = ", i, "/", n_area)
        if (length(areaList[[i]]) > 1) {
            count_new[, i] <- rowSums(seu[["RNA"]]@counts[, areaList[[i]]])
        } else {
            count_new[, i] <- seu[["RNA"]]@counts[, areaList[[i]]]
        }
        count_ct <- table(seu$merged_clusters[areaList[[i]]])
        domains <- names(which(count_ct == max(count_ct)))
        merged_clusters[i] <- domains[sample(seq_along(domains), 1)]
        feature_img_new[i, ] <- colMeans(feature_img[areaList[[i]], , drop = FALSE])
    }
    rm(seu)
    meta_data <- data.frame(
        row = pos_new[, 1], col = pos_new[, 2],
        merged_clusters = merged_clusters
    )
    row.names(meta_data) <- colnames(count_new)
    seu_new <- CreateSeuratObject(counts = as.sparse(count_new), meta.data = meta_data)
    list(
        seu = seu_new,
        feature_img = feature_img_new
    )
}

GSEA_enrich <- function(seu, pathways) {
    library(GSVA)
    a <- GSVA::gsva(
        seu@assays$RNA@data, pathways, method = "ssgsea", kcdf = "Gaussian",
        verbose = TRUE
    )

    pathway_gsea <- as.matrix(t(a))

    return(pathway_gsea)
}

GSVA_enrich <- function(seu, pathways) {
    library(GSVA)
    a <- GSVA::gsva(
        seu@assays$RNA@data, pathways, method = "gsva",
        verbose = TRUE
    )

    pathway_gsva <- as.matrix(t(a))

    return(pathway_gsva)
}

AUCell_enrich <- function(seu, pathways) {
    library(AUCell)

    pathway_AUCell <- AUCell_run(seu@assays$RNA@data, pathways)

    return(as.matrix(t(pathway_AUCell@assays@data$AUC)))
}

VAM_enrich <- function(seu, pathways) {
    library(VAM)

    gene.set.collection <- createGeneSetCollection(
        gene.ids = rownames(seu), gene.set.collection = pathways
    )
    obj <- vamForSeurat(
        seurat.data = seu,
        gene.set.collection = gene.set.collection,
        center = F, gamma = T, sample.cov = F, return.dist = F)
    vam.res <- as.matrix(t(obj@assays$VAMcdf@data))

    return(vam.res)
}













