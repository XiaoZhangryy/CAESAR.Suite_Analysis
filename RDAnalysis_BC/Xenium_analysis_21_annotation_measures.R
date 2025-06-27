rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

setwd(work_wd)

load("cols_cts_q.rda")
load("markerList.rda")

i <- as.integer(commandArgs(TRUE))

ana_wd <- paste0(work_wd, "BC", i, "/")
setwd(ana_wd)

annotation_wd <- paste0(ana_wd, "annotation_results/")
dir.create(annotation_wd, showWarnings = FALSE)

################################################################################
# 2. benchmark annotation
################################################################################
load("seu.rda")

load_results <- function(sample, annotation_wd = annotation_wd) {
    result_file_name0 <- paste0(
        annotation_wd, gsub(" ", "_", sample), "_CAESAR.rda"
    )
    load(result_file_name0)
    anno.df_CAESAR <- anno.df
    
    result_file_name0 <- paste0(
        annotation_wd, gsub(" ", "_", sample), "_CelliD.rda"
    )
    load(result_file_name0)
    anno.df_CelliD <- anno.df

    result_file_name1 <- paste0(
        annotation_wd, gsub(" ", "_", sample),
        "_Seurat_scmap_SingleR_scPred.rda"
    )
    load(result_file_name1)
    anno.df_benchmark <- anno.df[, paste0(c(
        "Seurat", "Seuratunasg", "Seuratconf",
        "scmap", "scmapunasg", "scmapconf",
        "SingleR", "SingleRunasg", "SingleRconf",
        "scPred", "scPredunasg", "scPredconf"
    ), "_", gsub(" ", "_", sample))]

    anno.df <- cbind(anno.df_CAESAR, anno.df_benchmark)
    anno.df <- cbind(anno.df, anno.df_CelliD)
    return(anno.df)
}

true_y <- seu$RCTD_first
cts_in_ref_all <- unique(unlist(lapply(markerList, names)))

## --------------------
## ACC
## --------------------

res_acc_pou <- lapply(names(markerList), function(sample, true_y, markerList, annotation_wd) {
    anno.df <- load_results(sample, annotation_wd)

    cts_in_ref <- names(markerList[[sample]])

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
    methods <- c(methods, paste0(methods, "unasg"))

    results <- sapply(methods, function(method) {
        mymeasures(anno.df[, paste0(method, "_", gsub(" ", "_", sample))],
        true_y, cts_in_ref)
    })
    colnames(results) <- methods
    return(results)
}, true_y = true_y, markerList = markerList, annotation_wd = annotation_wd)
names(res_acc_pou) <- names(markerList)


methods <- c("iCAESAR", "iCAESARunasg")
res_acc_pou_iCAESAR <- sapply(methods, function(method, true_y, cts_in_ref_all, annotation_wd) {
    result_file_name0 <- paste0(annotation_wd, "iCAESAR.rda")
    load(result_file_name0)

    mymeasures(anno.df[, method], true_y, cts_in_ref_all)
}, true_y = true_y, cts_in_ref_all = cts_in_ref_all, annotation_wd = annotation_wd)


## --------------------
## F1
## --------------------

res_Precision_Recall_F1 <- lapply(names(markerList), function(sample, true_y, markerList, annotation_wd) {
    anno.df <- load_results(sample, annotation_wd)

    cts_in_ref <- names(markerList[[sample]])

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")
    methods <- c(methods, paste0(methods, "unasg"))

    results <- lapply(methods, function(method) {
        calculate_f1_scores(anno.df[, paste0(method, "_", gsub(" ", "_", sample))],
        true_y, cts_in_ref)
    })
    names(results) <- methods
    return(results)
}, true_y = true_y, markerList = markerList, annotation_wd = annotation_wd)
names(res_Precision_Recall_F1) <- names(markerList)



methods <- c("iCAESAR", "iCAESARunasg")
res_Precision_Recall_F1_iCAESAR <- lapply(methods, function(method, true_y, cts_in_ref_all, annotation_wd) {
    result_file_name0 <- paste0(annotation_wd, "iCAESAR.rda")
    load(result_file_name0)

    calculate_f1_scores(anno.df[, method], true_y, cts_in_ref_all)
}, true_y = true_y, cts_in_ref_all = cts_in_ref_all, annotation_wd = annotation_wd)
names(res_Precision_Recall_F1_iCAESAR) <- methods



res_anno_measures <- list(
    acc_pou = res_acc_pou,
    acc_pou_iCAESAR = res_acc_pou_iCAESAR,
    f1score = res_Precision_Recall_F1,
    f1score_iCAESAR = res_Precision_Recall_F1_iCAESAR,
    n = table(true_y)
)

save(res_anno_measures, file = paste0(ana_wd, "res_anno_measures.rda"))


rm(
    res_acc_pou, res_acc_pou_iCAESAR, res_Precision_Recall_F1,
    res_Precision_Recall_F1_iCAESAR, res_anno_measures
)




## --------------------
## ASW and SigScore
## --------------------

res_asw <- c(
    CAESAR = asw(Embeddings(seu, "caesar"), seu$RCTD_first),
    CelliD = asw(Embeddings(seu, "mca"), seu$RCTD_first)
)
save(res_asw, file = "res_asw.rda")


load("degs_RCTD_first.rda")
topdegs <- degs %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

res_SigScore <- c(
    CAESAR = coembed_score(seu, "caesar", "RCTD_first", topdegs),
    CelliD = coembed_score(seu, "mca", "RCTD_first", topdegs)
)
save(res_SigScore, file = "res_SigScore.rda")


################################################################################
# 3. control pou 0.05
################################################################################
pou_th <- 0.05


## --------------------
## ACC
## --------------------

res_acc_pou <- lapply(names(markerList), function(
    sample, true_y, markerList, pou_th, annotation_wd) {
    anno.df <- load_results(sample, annotation_wd)

    cts_in_ref <- names(markerList[[sample]])

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")

    results <- sapply(methods, function(method) {
        pred_y <- pred_control_pou(
            anno.df,
            paste0(method, "unasg", "_", gsub(" ", "_", sample)),
            paste0(method, "_", gsub(" ", "_", sample)),
            paste0(method, "conf", "_", gsub(" ", "_", sample)),
            pou_th
        )

        mymeasures(pred_y, true_y, cts_in_ref)
    })
    colnames(results) <- methods
    return(results)
}, true_y = true_y, markerList = markerList,
pou_th = pou_th, annotation_wd = annotation_wd)
names(res_acc_pou) <- names(markerList)


methods <- c("iCAESAR")
res_acc_pou_iCAESAR <- sapply(methods, function(method, true_y, cts_in_ref_all, pou_th, annotation_wd) {
    result_file_name0 <- paste0(annotation_wd, "iCAESAR.rda")
    load(result_file_name0)

    pred_y <- pred_control_pou(
        anno.df,
        paste0(method, "unasg"),
        paste0(method),
        paste0(method, "conf"),
        pou_th
    )
    mymeasures(pred_y, true_y, cts_in_ref_all)
}, true_y = true_y, cts_in_ref_all = cts_in_ref_all,
pou_th = pou_th, annotation_wd = annotation_wd)



## --------------------
## F1
## --------------------

res_Precision_Recall_F1 <- lapply(names(markerList), function(
    sample, true_y, markerList, pou_th, annotation_wd) {
    anno.df <- load_results(sample, annotation_wd)

    cts_in_ref <- names(markerList[[sample]])

    methods <- 
        c("CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD")

    results <- lapply(methods, function(method) {
        pred_y <- pred_control_pou(
            anno.df,
            paste0(method, "unasg", "_", gsub(" ", "_", sample)),
            paste0(method, "_", gsub(" ", "_", sample)),
            paste0(method, "conf", "_", gsub(" ", "_", sample)),
            pou_th
        )

        calculate_f1_scores(pred_y, true_y, cts_in_ref)
    })
    names(results) <- methods
    return(results)
}, true_y = true_y, markerList = markerList,
pou_th = pou_th, annotation_wd = annotation_wd)
names(res_Precision_Recall_F1) <- names(markerList)


methods <- c("iCAESAR")
res_Precision_Recall_F1_iCAESAR <- lapply(methods, function(
    method, true_y, cts_in_ref_all, pou_th, annotation_wd) {
    result_file_name0 <- paste0(annotation_wd, "iCAESAR.rda")
    load(result_file_name0)

    pred_y <- pred_control_pou(
        anno.df,
        paste0(method, "unasg"),
        paste0(method),
        paste0(method, "conf"),
        pou_th
    )
    calculate_f1_scores(pred_y, true_y, cts_in_ref_all)
}, true_y = true_y, cts_in_ref_all = cts_in_ref_all,
pou_th = pou_th, annotation_wd = annotation_wd)
names(res_Precision_Recall_F1_iCAESAR) <- methods


res_anno_measures <- list(
    acc_pou = res_acc_pou,
    acc_pou_iCAESAR = res_acc_pou_iCAESAR,
    f1score = res_Precision_Recall_F1,
    f1score_iCAESAR = res_Precision_Recall_F1_iCAESAR,
    n = table(true_y)
)

save(res_anno_measures, file = paste0(ana_wd, "res_anno_measures_control_pou_0.05.rda"))




rm(
    res_acc_pou, res_acc_pou_iCAESAR, res_Precision_Recall_F1,
    res_Precision_Recall_F1_iCAESAR, res_anno_measures
)



