rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/MsHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

load("cols_cts_q.rda")


################################################################################
# 1. measures
################################################################################
methods <- c(
    "CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD"
)
custom_labels <- c(
    "CAESAR", "Seurat", "scmap", "SingleR", "scPred", "Cell-ID"
)
names(custom_labels) <- methods
custom_colors <- c(
    "#FF8C00", "#699ECA", "#2E8857", 
    "#9B59B6", "#008B8B", "#B8860B", "#556B8C"
)
names(custom_colors) <- methods


slice_colors <- c(
    "HCC1" = "#699ECA", "HCC2" = "#FF8C00",
    "HCC3" = "#F898CB", "HCC4" = "#4DAF4A"
)
slice_shapes <- c("HCC1" = 15, "HCC2" = 16, "HCC3" = 17, "HCC4" = 18)

#===============================================================================
# 1.1 plot measures - without account for unassigned
#===============================================================================
library(reshape2)
library(patchwork)
setwd(figs_wd)

plot_anno_HCC <- function(
    acc_melt, ylab, boxmethods, slice_colors = slice_colors, slice_shapes = slice_shapes) {
    plot1 <- ggplot(acc_melt, aes(x = Method, y = value)) +
        geom_boxplot(
            data = subset(acc_melt, Method %in% boxmethods),
            outlier.shape = NA, position = position_dodge(width = 0.95)
        ) + # Box plot with   adjusted position
        geom_jitter(position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.2), size = 1.5, aes(color = Slice, shape = Slice)) + # Scatter points with jitter
        scale_color_manual(values = slice_colors) + # Apply custom colors
        scale_shape_manual(values = slice_shapes) + # Apply custom shapes
        scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "white", color = "black", linewidth = 1), # Background with border
            panel.grid.major = element_blank(), # Remove major grid lines
            panel.grid.minor = element_blank(), # Remove minor grid lines
            axis.title.x = element_blank(), # Remove X-axis title
            axis.title.y = element_text(face = "bold", size = 20), # Y-axis title
            legend.position = "bottom",
            legend.text = element_text(face = "bold", size = 15),
            legend.title = element_text(face = "bold", size = 20, margin = margin(r = 10)),
            axis.text.y = element_text(face = "bold", size = 15),
            axis.text.x = element_text(face = "bold", size = 15)
        ) + # Position legend at the bottom
        ylim(0, 1) +
        labs(title = NULL, x = NULL, y = ylab)

    plot1
}


## --------------------
## ACC
## --------------------

acc_melt <- Reduce(rbind, lapply(1:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "HCC", i, "/")

    res_file <- paste0(ana_wd, "res_acc_pou_f1_7methods.rda")
    load(res_file)

    measures <- t(res_acc_pou["ACC", methods])
    dimnames(measures) <- list(
        Sample = "Visium",
        Method = methods
    )

    res <- reshape2::melt(measures)

    cbind(res, Slice = paste0("HCC", i))
}, methods = methods))
acc_melt$Method <- factor(acc_melt$Method, levels = methods)
acc_melt$Slice <- factor(acc_melt$Slice, levels = paste0("HCC", 1:4))

save(acc_melt, file = "acc_7methods_melt.rda")

plot0 <- plot_anno_HCC(acc_melt, "ACC", methods, slice_colors, slice_shapes)
ggsave(
    file = "ACC_7methods.png",
    plot = plot0, width = 9, height = 4, units = "in", dpi = 200
)





## --------------------
## F1
## --------------------
f1_melt <- Reduce(rbind, lapply(1:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "HCC", i, "/")

    res_file <- paste0(ana_wd, "res_acc_pou_f1_7methods.rda")
    load(res_file)

    measures <- t(f1score[methods])
    dimnames(measures) <- list(
        Sample = "Visium",
        Method = methods
    )

    res <- reshape2::melt(measures)

    cbind(res, Slice = paste0("HCC", i))
}, methods = methods))
f1_melt$Method <- factor(f1_melt$Method, levels = methods)
f1_melt$Slice <- factor(f1_melt$Slice, levels = paste0("HCC", 1:4))

save(f1_melt, file = "f1_7methods_melt.rda")

plot1 <- plot_anno_HCC(f1_melt, "F1", methods, slice_colors, slice_shapes)
ggsave(
    file = "F1_7methods.png",
    plot = plot1, width = 9, height = 4, units = "in", dpi = 200
)







#===============================================================================
# 1.2 control pou 0.05
#===============================================================================

# --------------------
### ACC
# --------------------

acc_melt <- Reduce(rbind, lapply(1:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "HCC", i, "/")

    res_file <- paste0(ana_wd, "res_acc_pou_f1_7methods_control_pou_0.05.rda")
    load(res_file)

    measures <- t(res_acc_pou["ACC", methods])
    dimnames(measures) <- list(
        Sample = "Visium",
        Method = methods
    )

    res <- reshape2::melt(measures)

    cbind(res, Slice = paste0("HCC", i))
}, methods = methods))
acc_melt$Method <- factor(acc_melt$Method, levels = methods)
acc_melt$Slice <- factor(acc_melt$Slice, levels = paste0("HCC", 1:4))

save(acc_melt, file = "acc_control_pou_0.05_7methods_melt.rda")


plot2 <- plot_anno_HCC(acc_melt, "ACC", methods, slice_colors, slice_shapes)
ggsave(
    file = "ACC_control_pou_0.05_7methods.png",
    plot = plot2, width = 9, height = 4, units = "in", dpi = 200
)




## --------------------
## F1
## --------------------
f1_melt <- Reduce(rbind, lapply(1:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "HCC", i, "/")

    res_file <- paste0(ana_wd, "res_acc_pou_f1_7methods_control_pou_0.05.rda")
    load(res_file)

    measures <- t(f1score[methods])
    dimnames(measures) <- list(
        Sample = "Visium",
        Method = methods
    )

    res <- reshape2::melt(measures)

    cbind(res, Slice = paste0("HCC", i))
}, methods = methods))
f1_melt$Method <- factor(f1_melt$Method, levels = methods)
f1_melt$Slice <- factor(f1_melt$Slice, levels = paste0("HCC", 1:4))

save(f1_melt, file = "f1_control_pou_0.05_7methods_melt.rda")

plot1 <- plot_anno_HCC(f1_melt, "F1", methods, slice_colors, slice_shapes)
ggsave(
    file = "F1_control_pou_0.05_7methods.png",
    plot = plot1, width = 9, height = 4, units = "in", dpi = 200
)




















































