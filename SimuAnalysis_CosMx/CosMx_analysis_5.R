################################################################################
# plot measures
################################################################################
rm(list = ls())
library(abind)
library(reshape2)
library(ggplot2)
library(dplyr)
work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMx/"
figs_wd <- paste0(work_wd, "Figs_all/")


sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)

enrich_melt <- Reduce(rbind, lapply(2:4, function(i) {
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(ana_wd, "celltypes.rda"))
    load(paste0(ana_wd, "enrich_measures_all.rda"))
    res_auc <- t(res_auc)
    dimnames(res_auc) <- list(
        cts = celltypes,
        Method = c("CAESAR", "Cell-ID", "GSDensity")
    )
    res <- melt(res_auc)
    cbind(res, Sample = sample_names[i])
}))
enrich_melt$Method <- factor(enrich_melt$Method, levels = c("CAESAR", "Cell-ID", "GSDensity"))
enrich_melt$Sample <- factor(enrich_melt$Sample, levels = sample_names[2:4])


save(enrich_melt, file = paste0(work_wd, "Figs_all/enrich_melt.rda"))

custom_labels <- c("CAESAR" = "CAESAR", "Cell-ID" = "Cell-ID", "GSDensity" = "GSDensity")
custom_colors <- c("CAESAR" = "#FF8C00", "Cell-ID" = "#699ECA", "GSDensity" = "#4DAF4A")
enrich_facet_plot <- ggplot(enrich_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1.0), limits = c(0.25, 1.0)) +
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
    # ylim(0.25, 1) +
    labs(title = NULL, x = NULL, y = "AUC Score")

ggsave(
    file = paste0(figs_wd, "Enrich_facet_plot.png"),
    plot = enrich_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)














