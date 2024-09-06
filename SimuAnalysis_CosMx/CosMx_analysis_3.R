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

sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)

################################################################################
# 2. summary and plot simulation results
################################################################################
sample_used <- c(
    "Lung5_Rep2", "Lung13", "Lung12"
)



## ACC
acc_melt <- Reduce(rbind, lapply(2:4, function(i) {
    # each query slice
    ana_wd <- paste0(work_wd, "Analysis_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "measures_fov", qfov, ".rda")
        load(res_file)
        res1 <- t(sapply(measures, function(dat) dat["ACC", c("CAESARunasg", "CelliDunasg")]))
        res1 <- cbind(iCAESAR_measures["ACC", 2], res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", "CAESAR", "CelliD"),
        qfov = fovs
    )
    res <- melt(measures_array)
    cbind(res, Sample = sample_names[i])
}))
acc_melt$Method <- factor(acc_melt$Method, levels = c("iCAESAR", "CAESAR", "CelliD"))
acc_melt$Sample <- factor(acc_melt$Sample, levels = sample_used)

save(acc_melt, file = paste0(work_wd, "Figs_all/acc_melt.rda"))

custom_labels <- c("iCAESAR" = "iCAESAR", "CAESAR" = "CAESAR", "CelliD" = "Cell-ID")
custom_colors <- c("iCAESAR" = "#C74546", "CAESAR" = "#FF8C00", "CelliD" = "#699ECA")
acc_facet_plot <- ggplot(acc_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(values = custom_colors) + # Apply custom colors
    scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
        legend.position = "right",
        axis.title.y = element_text(face = "bold", size = 12), # Y-axis title
        axis.text.y = element_text(face = "bold", size = 10),
        # axis.ticks.y = element_blank(), # Remove y-axis ticks for all but leftmost plot
        strip.background = element_blank(), # Remove background of facet labels
        strip.text = element_text(face = "bold"), # Bold facet labels
        panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
        plot.margin = unit(c(1, 1, 1, 1), "cm"), # Adjust plot margins
        panel.grid.major.y = element_blank(), # Remove y-axis grid lines
        panel.grid.minor.y = element_blank(), # Remove y-axis minor grid lines
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.title.x = element_blank(), # Remove X-axis title
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 10) # Adjust size of facet labels
    ) +
    ylim(0, 1) +
    labs(title = NULL, x = NULL, y = "ACC")

ggsave(
    file = paste0(figs_wd, "ACC_facet_plot.png"),
    plot = acc_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)

## ASW
asw_melt <- Reduce(rbind, lapply(2:4, function(i) {
    # each query slice
    ana_wd <- paste0(work_wd, "Analysis_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "embedding_measures_fov", qfov, ".rda")
        load(res_file)
        res_asw
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 2)
    dimnames(measures_array) <- list(
        Method = c("CAESAR", "CelliD"),
        qfov = fovs
    )
    res <- melt(measures_array)
    cbind(res, Sample = sample_names[i])
}))
asw_melt$Method <- factor(asw_melt$Method, levels = c("CAESAR", "CelliD"))
asw_melt$Sample <- factor(asw_melt$Sample, levels = sample_used)

save(asw_melt, file = paste0(work_wd, "Figs_all/asw_melt.rda"))


custom_labels <- c("CAESAR" = "CAESAR", "CelliD" = "Cell-ID")
custom_colors <- c("CAESAR" = "#FF8C00", "CelliD" = "#699ECA")
asw_facet_plot <- ggplot(asw_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(values = custom_colors) + # Apply custom colors
    scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
        legend.position = "right",
        axis.title.y = element_text(face = "bold", size = 12), # Y-axis title
        axis.text.y = element_text(face = "bold", size = 10),
        # axis.ticks.y = element_blank(), # Remove y-axis ticks for all but leftmost plot
        strip.background = element_blank(), # Remove background of facet labels
        strip.text = element_text(face = "bold"), # Bold facet labels
        panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
        plot.margin = unit(c(1, 1, 1, 1), "cm"), # Adjust plot margins
        panel.grid.major.y = element_blank(), # Remove y-axis grid lines
        panel.grid.minor.y = element_blank(), # Remove y-axis minor grid lines
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.title.x = element_blank(), # Remove X-axis title
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 10) # Adjust size of facet labels
    ) +
    ylim(-0.25, 0.5) +
    labs(title = NULL, x = NULL, y = "ASW")


ggsave(
    file = paste0(figs_wd, "ASW_facet_plot.png"),
    plot = asw_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


## signature score
sigscore_melt <- Reduce(rbind, lapply(2:4, function(i) {
    # each query slice
    ana_wd <- paste0(work_wd, "Analysis_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "coembedding_measures_fov", qfov, ".rda")
        load(res_file)
        res_coembed_score
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 2)
    dimnames(measures_array) <- list(
        Method = c("CAESAR", "CelliD"),
        qfov = fovs
    )
    res <- melt(measures_array)
    cbind(res, Sample = sample_names[i])
}))
sigscore_melt$Method <- factor(sigscore_melt$Method, levels = c("CAESAR", "CelliD"))
sigscore_melt$Sample <- factor(sigscore_melt$Sample, levels = sample_used)

save(sigscore_melt, file = paste0(work_wd, "Figs_all/sigscore_melt.rda"))

custom_labels <- c("CAESAR" = "CAESAR", "CelliD" = "Cell-ID")
custom_colors <- c("CAESAR" = "#FF8C00", "CelliD" = "#699ECA")
sigscore_facet_plot <- ggplot(sigscore_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(values = custom_colors) + # Apply custom colors
    scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1.0), limits = c(0.25, 1.0)) +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
        legend.position = "right",
        axis.title.y = element_text(face = "bold", size = 12), # Y-axis title
        axis.text.y = element_text(face = "bold", size = 10),
        # axis.ticks.y = element_blank(), # Remove y-axis ticks for all but leftmost plot
        strip.background = element_blank(), # Remove background of facet labels
        strip.text = element_text(face = "bold"), # Bold facet labels
        panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
        plot.margin = unit(c(1, 1, 1, 1), "cm"), # Adjust plot margins
        panel.grid.major.y = element_blank(), # Remove y-axis grid lines
        panel.grid.minor.y = element_blank(), # Remove y-axis minor grid lines
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.title.x = element_blank(), # Remove X-axis title
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 10) # Adjust size of facet labels
    ) +
    # ylim(0.25, 1) +
    labs(title = NULL, x = NULL, y = "Signature Score")


ggsave(
    file = paste0(figs_wd, "sigscore_facet_plot_plot.png"),
    plot = sigscore_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


















