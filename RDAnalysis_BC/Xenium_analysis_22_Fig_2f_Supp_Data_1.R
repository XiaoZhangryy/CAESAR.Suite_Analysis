rm(list = ls())
################################################################################
# Figure 2.f
################################################################################

work_wd <- "/share/analysisdata/zhangx/CAESAR/XeniumRDA/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")

################################################################################
# 1. plot measures
################################################################################
methods <- c(
    "CAESAR", "Seurat", "scmap", "SingleR", "scPred", "CelliD"
)
custom_labels <- c(
    "iCAESAR", "CAESAR", "Seurat", "scmap", "SingleR", "scPred", "Cell-ID"
)
names(custom_labels) <- c("iCAESAR", methods)
custom_colors <- c(
    "#C74546", "#FF8C00", "#699ECA", "#2E8857", 
    "#9B59B6", "#008B8B", "#B8860B", "#556B8C"
)
names(custom_colors) <- c("iCAESAR", methods)


slice_colors <- c(
    "BC1" = "#699ECA", "BC2" = "#FF8C00",
    "BC3" = "#F898CB", "BC4" = "#4DAF4A"
)
slice_shapes <- c("BC1" = 15, "BC2" = 16, "BC3" = 17, "BC4" = 18)



plot_methods <- c(
  "iCAESAR", "CAESAR", "Seurat", "scmap", "scPred", "SingleR", "CelliD"
)
custom_labels <- custom_labels[plot_methods]
custom_colors <- custom_colors[plot_methods]


plot_anno_Xenium <- function(
    acc_melt, ylab, boxmethods, slice_colors = slice_colors, slice_shapes = slice_shapes) {
    plot1 <- ggplot(acc_melt, aes(x = Method, y = value)) +
        geom_boxplot(
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








#===============================================================================
# 1.1 plot measures - fig. 2.f
#===============================================================================
library(reshape2)
library(patchwork)
setwd(figs_wd)

load("acc_control_pou_0.05_7methods_melt.rda")

acc_melt <- acc_melt[acc_melt$Method %in% plot_methods, ]
acc_melt$Method <- factor(acc_melt$Method, levels = plot_methods)

plot1 <- plot_anno_Xenium(acc_melt, "ACC", methods, slice_colors, slice_shapes) +
  theme(
    axis.title.y = element_text(face = "bold", size = 20), # Y-axis title
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.text.x = element_text(face = "bold", size = 15, angle = 30, hjust = 1)
  )
ggsave(
  file = "ACC_7methods_combineBC.png",
  plot = plot1, width = 6, height = 4, units = "in", dpi = 200
)
ggsave(
  file = "ACC_7methods_combineBC_nolegend.png",
  plot = plot1 + theme(legend.position = "none"), width = 6, height = 3.5, units = "in", dpi = 200
)



## --------------------
## 1.2 plot embedding
## --------------------
res_asw <- sapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_asw.rda"))
    res_asw
})
colnames(res_asw) <- paste0("BC", 1:4)
rownames(res_asw) <- c("CAESAR", "CelliD")

res_SigScore <- sapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_SigScore.rda"))
    res_SigScore
})
colnames(res_SigScore) <- paste0("BC", 1:4)
rownames(res_SigScore) <- c("CAESAR", "CelliD")

save(res_asw, res_SigScore, file = "ASW_SigScore.rda")




plot_embed_Xenium <- function(embed, ylab) {
    data_df <- as.data.frame(embed)
    data_df$Slice <- rownames(data_df)

    data_melted <- melt(data_df, id.vars = "Slice", variable.name = "Method", value.name = "Result")

    custom_colors <- c(
        "BC1" = "#699ECA", "BC2" = "#FF8C00",
        "BC3" = "#F898CB", "BC4" = "#4DAF4A"
    )
    custom_shapes <- c("BC1" = 15, "BC2" = 16, "BC3" = 17, "BC4" = 18)
    custom_labels <- c("CAESAR" = "CAESAR", "CelliD" = "Cell-ID")

    plot1 <- ggplot(data_melted, aes(x = Method, y = Result)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) + # Box plot with adjusted position
        geom_jitter(position = position_jitter(width = 0.2), size = 3, aes(color = Slice, shape = Slice)) + # Scatter points with jitter
        scale_color_manual(values = custom_colors) + # Apply custom colors
        scale_shape_manual(values = custom_shapes) + # Apply custom shapes
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
        labs(title = NULL, x = NULL, y = ylab)

    return(plot1)
}

plot2 <- plot_embed_Xenium(t(res_asw), "ASW")
plot3 <- plot_embed_Xenium(t(res_SigScore), "SigScore")

library(patchwork)
combined_plot <- plot2 + plot_spacer() + plot3 + plot_layout(widths = c(10, 0.5, 10), guides = "collect") & theme(legend.position = "bottom")

# combined_plot <- plot1 + plot_spacer() + plot2 + plot_layout(widths = c(10, 0.5, 10), guides = "collect") & theme(legend.position = "none")

ggsave(
    file = "ASW_SigScore.png",
    plot = combined_plot, width = 6, height = 3, units = "in", dpi = 200
)










library(patchwork)
plot1 <- plot1 + theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 15, angle = 25, hjust = 1)
)


combined_plot2 <- plot1 / combined_plot

ggsave(
    file = "ACC_ASW_SigScore.png",
    plot = combined_plot2, width = 6, height = 6, units = "in", dpi = 200
)







#===============================================================================
# 1.2 plot measures with unassigned
#===============================================================================
acc_melt <- Reduce(rbind, lapply(1:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "BC", i, "/")

    res_file <- paste0(ana_wd, "res_anno_measures.rda")
    load(res_file)

    measures <- t(sapply(res_anno_measures$acc_pou, function(dat) dat["ACC", paste0(methods, "unasg")]))
    dimnames(measures) <- list(
        Sample = rownames(measures),
        Method = methods
    )

    res <- reshape2::melt(measures)

    res <- rbind(
        data.frame(
            "Sample" = "Sample 0",
            "Method" = "iCAESAR",
            "value" = res_anno_measures$acc_pou_iCAESAR["ACC", "iCAESARunasg"]
        ), res
    )

    cbind(res, Slice = paste0("BC", i))
}, methods = methods))
acc_melt$Method <- factor(acc_melt$Method, levels = c("iCAESAR", methods))
acc_melt$Slice <- factor(acc_melt$Slice, levels = paste0("BC", 1:4))

save(acc_melt, file = "acc_unasg_7methods_melt.rda")


load("acc_unasg_7methods_melt.rda")

plot_anno_Xenium <- function(
    acc_melt, ylab, boxmethods, slice_colors = slice_colors, slice_shapes = slice_shapes) {
    plot1 <- ggplot(acc_melt, aes(x = Method, y = value, color = Slice, shape = Slice)) +
        geom_boxplot(
            data = subset(acc_melt, Method %in% boxmethods),
            outlier.shape = NA, position = position_dodge(width = 0.95)
        ) + # Box plot with   adjusted position
        geom_jitter(position = position_jitterdodge(dodge.width = 0.95, jitter.width = 0.2), size = 1.5) + # Scatter points with jitter
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


acc_melt <- acc_melt[acc_melt$Method %in% plot_methods, ]
acc_melt$Method <- factor(acc_melt$Method, levels = plot_methods)


plot0 <- plot_anno_Xenium(acc_melt, "ACC", methods, slice_colors, slice_shapes)
ggsave(
    file = "ACC_unasg_7methods.png",
    plot = plot0, width = 9, height = 5, units = "in", dpi = 200
)


################################################################################
# 2. markerList to csv
################################################################################
load("markerList.rda")

markermatList <- lapply(markerList, function(marker) {
    sapply(marker, identity)
})

library(openxlsx)

wb <- createWorkbook()

for (i in seq_along(markermatList)) {
  # Add a sheet with the matrix name
  addWorksheet(wb, sheetName = names(markermatList)[i])
  
  # Write the matrix to the sheet
  writeData(wb, sheet = names(markermatList)[i], markermatList[[i]])
}

saveWorkbook(wb, file = "BC_markers.xlsx", overwrite = TRUE)

