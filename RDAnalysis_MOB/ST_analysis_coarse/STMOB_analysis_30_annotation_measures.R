rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/STMOB_coarse/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "STMOB/")
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results/")

load("seu.rda")

################################################################################
# 2. measures
################################################################################
manual_y <- function(manual_annotation) {
    manual_annotation <- as.character(manual_annotation)
    manual_annotation[manual_annotation == "GCL"] <- "GC"
    manual_annotation[manual_annotation == "MCL"] <- "M/TC"
    manual_annotation[manual_annotation == "ONL"] <- "OSNs"
    manual_annotation[manual_annotation == "GL"] <- "PGC"
    return(manual_annotation)
}

load_results <- function(annotation_wd = annotation_wd) {
    result_file_name1 <- paste0(
        annotation_wd, "CAESAR_Seurat_scmap_SingleR_scPred_CelliD.rda"
    )
    load(result_file_name1)
    return(anno.df)
}

anno.df <- load_results(annotation_wd)
true_y <- manual_y(seu$manual_annotation)


## --------------------
## ACC and F1
## --------------------
methods <- 
    c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
methods <- c(methods, paste0(methods, "unasg"))

res_acc_pou <- sapply(methods, function(method, anno.df, true_y) {
    pred_y <- anno.df[, method]
    c(
        "ACC" = mean(true_y == pred_y),
        "POU" = mean(pred_y == "unassigned")
    )
}, anno.df = anno.df, true_y = true_y)

f1score_list <- lapply(methods, function(method, anno.df, true_y) {
    pred_y <- anno.df[, method]
    calculate_f1_scores(pred_y, true_y, c("GC", "M/TC", "OSNs", "PGC"))
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
    calculate_f1_scores(pred_y, true_y, c("GC", "M/TC", "OSNs", "PGC"))
}, anno.df = anno.df, true_y = true_y, pou_th = pou_th)

f1score <- sapply(f1score_list, function(F1res, true_y) {
    n_cts <- table(true_y)
    nn <- n_cts[colnames(F1res)]
    nn <- nn / sum(nn)
    F1res["F1", ] %*% nn
}, true_y = true_y)
names(f1score) <- methods


save(res_acc_pou, f1score_list, f1score, file = "res_acc_pou_f1_7methods_control_pou_0.05.rda")




