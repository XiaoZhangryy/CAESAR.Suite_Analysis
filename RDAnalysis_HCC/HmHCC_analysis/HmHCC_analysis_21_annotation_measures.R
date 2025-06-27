rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/HmHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

load("cols_cts_q.rda")

i <- as.integer(commandArgs(TRUE))

ana_wd <- paste0(work_wd, "HCC", i, "/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results/")

load("seu.rda")

################################################################################
# 1. measures
################################################################################
load_results <- function(annotation_wd = annotation_wd) {
    result_file_name1 <- paste0(
        annotation_wd, "CAESAR_Seurat_scmap_SingleR_scPred_CelliD.rda"
    )
    load(result_file_name1)
    return(anno.df)
}

anno.df <- load_results(annotation_wd)
true_y <- seu$manual_labels

labels201 <- function(y) {
    y0 <- rep(1, length(y))
    y0[y == "Malignant cell"] <- 0
    y0[y == "HPC-like"] <- 0
    y0[y == "unassigned"] <- 0.5
    y0
}

## --------------------
## ACC and F1
## --------------------
methods <- 
    c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
methods <- c(methods, paste0(methods, "unasg"))

res_acc_pou <- sapply(methods, function(method, anno.df, true_y) {
    pred_y <- labels201(anno.df[, method])
    c(
        "ACC" = mean(true_y == pred_y),
        "POU" = mean(pred_y == "unassigned")
    )
}, anno.df = anno.df, true_y = true_y)

f1score_list <- lapply(methods, function(method, anno.df, true_y) {
    pred_y <- labels201(anno.df[, method])
    calculate_f1_scores(pred_y, true_y, unique(true_y))
}, anno.df = anno.df, true_y = true_y)

f1score <- sapply(f1score_list, function(F1res, true_y) {
    n_cts <- table(true_y)
    nn <- n_cts[colnames(F1res)]
    nn <- nn / sum(nn)
    F1res["F1", ] %*% nn
}, true_y = true_y)
names(f1score) <- methods


save(res_acc_pou, f1score_list, f1score, file = "res_acc_pou_f1_7methods.rda")



## --------------------
## control pou 0.05
## --------------------
pou_th <- 0.05


methods <- 
    c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")

res_acc_pou <- sapply(methods, function(method, anno.df, true_y, pou_th) {
    pred_y <- pred_control_pou(
        anno.df,
        paste0(method, "unasg"),
        paste0(method),
        paste0(method, "conf"),
        pou_th
    )
    pred_y <- labels201(pred_y)
    c(
        "ACC" = mean(true_y == pred_y),
        "POU" = mean(pred_y == "unassigned")
    )
}, anno.df = anno.df, true_y = true_y, pou_th = pou_th)

f1score_list <- lapply(methods, function(method, anno.df, true_y, pou_th) {
    pred_y <- pred_control_pou(
        anno.df,
        paste0(method, "unasg"),
        paste0(method),
        paste0(method, "conf"),
        pou_th
    )
    pred_y <- labels201(pred_y)
    calculate_f1_scores(pred_y, true_y, unique(true_y))
}, anno.df = anno.df, true_y = true_y, pou_th = pou_th)

f1score <- sapply(f1score_list, function(F1res, true_y) {
    n_cts <- table(true_y)
    nn <- n_cts[colnames(F1res)]
    nn <- nn / sum(nn)
    F1res["F1", ] %*% nn
}, true_y = true_y)
names(f1score) <- methods


save(res_acc_pou, f1score_list, f1score, file = "res_acc_pou_f1_7methods_control_pou_0.05.rda")








