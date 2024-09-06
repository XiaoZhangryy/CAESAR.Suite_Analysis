rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/Xenium/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")

################################################################################
# 1. measures
################################################################################
## calculate pou
calculate_pou <- function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    load("markerList.rda")
    load("cols_cts_q.rda")

    res_pou <- sapply(names(markerList), function(sample) {
        marker <- markerList[[sample]]

        # CAESAR
        methods <- paste0(c("CAESARunasg", "CelliDunasg"), "_", gsub(" ", "_", sample))
        results <- sapply(methods, function(method) {
            mean(as.character(seu@meta.data[, method]) == "unassigned")
        })
        names(results) <- c("CAESARunasg", "CelliDunasg")
        results
    })

    res_pou_iCAESAR <- mean(as.character(seu$iCAESARunasg) == "unassigned")

    save(res_pou, res_pou_iCAESAR, file = paste0(ana_wd, "res_pou.rda"))
}
sapply(1:4, calculate_pou)

recalculate_acc <- function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "seu.rda"))
    load("markerList.rda")
    load("cols_cts_q.rda")

    res_acc <- sapply(names(markerList), function(sample) {
        marker <- markerList[[sample]]

        cts_notin_sample <- setdiff(celltypes, names(marker))
        # SAX
        methods <- paste0(c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg"), "_", gsub(" ", "_", sample))
        results <- sapply(methods, function(method) {
            calculate_acc(seu@meta.data[, method], seu$RCTD_first, cts_notin_sample)
        })
        names(results) <- c("CAESAR", "CAESARunasg", "CelliD", "CelliDunasg")
        results
    })

    res_acc_iCAESAR <- sapply(c("iCAESAR", "iCAESARunasg"), function(method) {
        calculate_acc(seu@meta.data[, method], seu$RCTD_first, NULL)
    })
    names(res_acc_iCAESAR) <- c("iCAESAR", "iCAESARunasg")

    save(res_acc, res_acc_iCAESAR, file = paste0(ana_wd, "res_acc.rda"))
}
sapply(1:4, recalculate_acc)

## summary measures
res_pou <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_pou.rda"))
    t(res_pou)
})
names(res_pou) <- paste0("BC", 1:4)

res_ipou <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_pou.rda"))
    res_pou_iCAESAR
})
names(res_ipou) <- paste0("BC", 1:4)

res_acc <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_acc.rda"))
    t(res_acc)
})
names(res_acc) <- paste0("BC", 1:4)

res_iacc <- lapply(1:4, function(i) {
    ana_wd <- paste0(work_wd, "BC", i, "/")
    load(paste0(ana_wd, "res_acc.rda"))
    res_acc_iCAESAR
})
names(res_iacc) <- paste0("BC", 1:4)

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

save(res_pou, res_ipou, res_acc, res_iacc, res_asw, res_SigScore, file = "res_measures.rda")


#===============================================================================
# 1.1 plot measures - fig. 2.f
#===============================================================================
library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)

plot_anno_Xenium <- function(ACC, ACC_iSAX, ylab) {
    data_list <- lapply(names(ACC), function(slice) {
        df <- as.data.frame(ACC[[slice]])
        df$Slice <- slice
        df$Sample <- 1:nrow(df)
        df
    })

    data_combined <- do.call(rbind, data_list)
    data_melted <- melt(data_combined, id.vars = c("Slice", "Sample"), variable.name = "Method", value.name = "Result")

    if (!is.null(ACC_iSAX)) {
        data_melted_iSAX <- melt(do.call(rbind, ACC_iSAX), id.vars = c("Slice", "Sample"), variable.name = "Method", value.name = "Result")

        data_melted <- rbind(data_melted_iSAX, data_melted)
    }
    data_melted <- data_melted[data_melted$Method %in% c("iCAESARua", "CAESARua", "CelliDua"), ]
    data_melted$Method <- factor(data_melted$Method, levels = c("iCAESARua", "CAESARua", "CelliDua"))

    custom_colors <- c(
        "BC1" = "#699ECA", "BC2" = "#FF8C00",
        "BC3" = "#F898CB", "BC4" = "#4DAF4A"
    )
    custom_shapes <- c("BC1" = 15, "BC2" = 16, "BC3" = 17, "BC4" = 18)
    custom_labels <- c("iCAESARua" = "iCAESAR", "CAESARua" = "CAESAR", "CelliDua" = "Cell-ID")

    plot1 <- ggplot(data_melted, aes(x = Method, y = Result, color = Slice, shape = Slice)) +
        geom_boxplot(
            data = subset(data_melted, Method %in% c("CAESARua", "CelliDua")),
            outlier.shape = NA, position = position_dodge(width = 0.9)
        ) + # Box plot with   adjusted position
        geom_jitter(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2), size = 2) + # Scatter points with jitter
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
        ylim(0, 1) +
        labs(title = NULL, x = NULL, y = ylab)

    plot1
}

plot0 <- plot_anno_Xenium(res_acc, lapply(names(res_iacc), function(situation) {
  acc <- res_iacc[[situation]]
  df <- data.frame(CAESAR = acc[1], CAESARua = acc[2], Slice = situation, Sample = 0)
  colnames(df)[1:2] <- names(acc)
  df
}), "ACC")
ggsave(
    file = "ACC.png",
    plot = plot0, width = 6, height = 4, units = "in", dpi = 200
)

ggsave(
    file = "ACC_nolegend.png",
    plot = plot0 + theme(legend.position = "none"), width = 6, height = 3, units = "in", dpi = 200
)

## 1.2 plot embedding
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


