rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMxSimu/"
setwd(work_wd)
source("functions.R")

number <- commandArgs(TRUE) %>%
    as.integer() %>%
    assign4CosMx()
i <- number[1]
fov <- number[2]

figs_wd <- paste0(work_wd, "Figs/")
ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/fov", fov, "/")
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results_scenario3/")

seu_name_scenario3 <- paste0(ana_wd, "seu_fov", fov, "_scenario3.rda")
load(seu_name_scenario3)

## reference data
ref_sample <- "Lung5_Rep1"
ref_i <- 1
ref_process_wd <- paste0(work_wd, "processed_data_CosMx", ref_i, "/")
load(paste0(ref_process_wd, "fovs.rda"))
load(paste0(ref_process_wd, "q_est.rda"))
load(paste0(ref_process_wd, "cols.rda"))
cols0 <- cols

################################################################################
# 1. measures
################################################################################
cts_refList <- lapply(fovs, function(rfov) {
    load(paste0(ref_process_wd, "fov", rfov, "/sig_list_scRNA_fov", rfov, ".rda"))
    names(sig_list[[1]])
})
names(cts_refList) <- fovs

################################################################################
# 2. without account for unassigned
################################################################################
load_results <- function(rfov, ref_i = ref_i, annotation_wd = annotation_wd) {
    result_file_name1 <- paste0(
        annotation_wd, "ref", ref_i, "fov", rfov,
        "_CAESAR_Seurat_scmap_SingleR_scPred_CelliD.rda")
    load(result_file_name1)
    return(anno.df)
}

true_y <- seu$merged_clusters

res_acc_pou <- lapply(fovs, function(rfov, true_y, cts_refList, ref_i, annotation_wd) {
    anno.df <- load_results(rfov, ref_i, annotation_wd)

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
    methods <- c(methods, paste0(methods, "unasg"))

    results <- sapply(methods, function(method) {
        mymeasures(anno.df[, paste0(method, "_ref", ref_i, "fov", rfov)],
        true_y, cts_refList[[rfov]])
    })
    colnames(results) <- methods
    return(results)
}, true_y = true_y, cts_refList = cts_refList, ref_i = ref_i,
annotation_wd = annotation_wd)
names(res_acc_pou) <- fovs

methods <- paste0(
    c("iCAESAR", "iCAESARunasg"), "_ref", ref_i
)
res_acc_pou_iCAESAR <- sapply(methods, function(method, cts_refList) {
    mymeasures(seu@meta.data[, method], seu$merged_clusters, unique(unlist(cts_refList)))
}, cts_refList = cts_refList)




## F1
res_Precision_Recall_F1 <- lapply(fovs, function(rfov, true_y, cts_refList, ref_i, annotation_wd) {
    anno.df <- load_results(rfov, ref_i, annotation_wd)

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
    methods <- c(methods, paste0(methods, "unasg"))

    results <- lapply(methods, function(method) {
        calculate_f1_scores(anno.df[, paste0(method, "_ref", ref_i, "fov", rfov)],
        true_y, cts_refList[[rfov]])
    })
    names(results) <- methods
    return(results)
}, true_y = true_y, cts_refList = cts_refList, ref_i = ref_i,
annotation_wd = annotation_wd)
names(res_Precision_Recall_F1) <- fovs


methods <- paste0(
    c("iCAESAR", "iCAESARunasg"), "_ref", ref_i
)
res_Precision_Recall_F1_iCAESAR <- lapply(methods, function(method, cts_refList) {
    calculate_f1_scores(seu@meta.data[, method], seu$merged_clusters, unique(unlist(cts_refList)))
}, cts_refList = cts_refList)
names(res_Precision_Recall_F1_iCAESAR) <- c("iCAESAR", "iCAESARunasg")



res_anno_measures_scenario3 <- list(
    acc_pou = res_acc_pou,
    acc_pou_iCAESAR = res_acc_pou_iCAESAR,
    f1score = res_Precision_Recall_F1,
    f1score_iCAESAR = res_Precision_Recall_F1_iCAESAR,
    n = table(seu$merged_clusters)
)

save(res_anno_measures_scenario3, file = paste0(ana_wd, "res_anno_measures_fov", fov, "_scenario3.rda"))





################################################################################
# 3. control pou 0.05
################################################################################
pou_th <- 0.05

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

res_acc_pou <- lapply(fovs, function(rfov, true_y, cts_refList, ref_i, annotation_wd, pou_th) {
    anno.df <- load_results(rfov, ref_i, annotation_wd)

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
    results <- sapply(methods, function(method) {
        pred_y <- pred_control_pou(
            anno.df,
            paste0(method, "unasg_ref", ref_i, "fov", rfov),
            paste0(method, "_ref", ref_i, "fov", rfov),
            paste0(method, "conf_ref", ref_i, "fov", rfov),
            pou_th
        )

        mymeasures(pred_y, true_y, cts_refList[[rfov]])
    })
    colnames(results) <- methods
    return(results)
}, true_y = true_y, cts_refList = cts_refList, ref_i = ref_i,
annotation_wd = annotation_wd, pou_th = pou_th)
names(res_acc_pou) <- fovs

methods <- c("iCAESAR")
res_acc_pou_iCAESAR <- sapply(methods, function(method, cts_refList, true_y, pou_th) {
    pred_y <- pred_control_pou(
        seu@meta.data,
        paste0(method, "unasg_ref", ref_i),
        paste0(method, "_ref", ref_i),
        paste0(method, "conf_ref", ref_i),
        pou_th
    )

    mymeasures(pred_y, true_y, unique(unlist(cts_refList)))
}, cts_refList = cts_refList, true_y = true_y, pou_th = pou_th)




## F1
res_Precision_Recall_F1 <- lapply(fovs, function(rfov, true_y, cts_refList, ref_i, annotation_wd, pou_th) {
    anno.df <- load_results(rfov, ref_i, annotation_wd)

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")

    results <- lapply(methods, function(method) {
        pred_y <- pred_control_pou(
            anno.df,
            paste0(method, "unasg_ref", ref_i, "fov", rfov),
            paste0(method, "_ref", ref_i, "fov", rfov),
            paste0(method, "conf_ref", ref_i, "fov", rfov),
            pou_th
        )

        calculate_f1_scores(pred_y, true_y, cts_refList[[rfov]])
    })
    names(results) <- methods
    return(results)
}, true_y = true_y, cts_refList = cts_refList, ref_i = ref_i,
annotation_wd = annotation_wd, pou_th = pou_th)
names(res_Precision_Recall_F1) <- fovs



methods <- c("iCAESAR")
res_Precision_Recall_F1_iCAESAR <- lapply(methods, function(method, cts_refList, true_y, pou_th) {
    pred_y <- pred_control_pou(
        seu@meta.data,
        paste0(method, "unasg_ref", ref_i),
        paste0(method, "_ref", ref_i),
        paste0(method, "conf_ref", ref_i),
        pou_th
    )

    calculate_f1_scores(pred_y, true_y, unique(unlist(cts_refList)))
}, cts_refList = cts_refList, true_y = true_y, pou_th = pou_th)
names(res_Precision_Recall_F1_iCAESAR) <- methods


res_anno_measures_scenario3 <- list(
    acc_pou = res_acc_pou,
    acc_pou_iCAESAR = res_acc_pou_iCAESAR,
    f1score = res_Precision_Recall_F1,
    f1score_iCAESAR = res_Precision_Recall_F1_iCAESAR,
    n = table(true_y)
)

save(res_anno_measures_scenario3, file = paste0(ana_wd, "res_anno_measures_control_pou_0.05_fov", fov, "_scenario3.rda"))









