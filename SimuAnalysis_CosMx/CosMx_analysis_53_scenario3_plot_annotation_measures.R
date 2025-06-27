rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/CosMxSimu/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

sample_names <- c(
    "Lung5_Rep1", "Lung5_Rep2", "Lung12", "Lung13"
)

################################################################################
# 2. summary and plot simulation results
################################################################################
sample_used <- c(
    "Lung5_Rep2", "Lung13", "Lung12"
)

# ==============================================================================
## without account for unassigned
# ==============================================================================
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

# --------------------
### ACC
# --------------------

acc_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_fov", qfov, "_scenario3.rda")
        load(res_file)
        res1 <- t(sapply(res_anno_measures_scenario3$acc_pou, function(dat) dat["ACC", methods]))
        res1 <- cbind(res_anno_measures_scenario3$acc_pou_iCAESAR["ACC", "iCAESAR_ref1"], res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
acc_melt$Method <- factor(acc_melt$Method, levels = c("iCAESAR", methods))
acc_melt$Sample <- factor(acc_melt$Sample, levels = sample_used)

save(acc_melt, file = paste0(figs_wd, "scenario3_acc_7methods_melt.rda"))



acc_facet_plot <- ggplot(acc_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    file = paste0(figs_wd, "scenario3_ACC_7methods_facet_plot.png"),
    plot = acc_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)

rm(acc_melt, acc_facet_plot)

# --------------------
### F1
# --------------------
f1_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_fov", qfov, "_scenario3.rda")
        load(res_file)

        n_cts <- res_anno_measures_scenario3$n

        res1 <- t(sapply(res_anno_measures_scenario3$f1score, function(dat) {
            sapply(dat[methods], function(F1res) {
                nn <- n_cts[colnames(F1res)]
                nn <- nn / sum(nn)
                F1res["F1", ] %*% nn
            })
        }))

        F1res <- res_anno_measures_scenario3$f1score_iCAESAR[["iCAESAR"]]
        nn <- n_cts[colnames(F1res)]
        nn <- nn / sum(nn)
        F1_iCAESAR <- as.numeric(F1res["F1", ] %*% nn)
        
        res1 <- cbind(F1_iCAESAR, res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
f1_melt$Method <- factor(f1_melt$Method, levels = c("iCAESAR", methods))
f1_melt$Sample <- factor(f1_melt$Sample, levels = sample_used)

save(f1_melt, file = paste0(figs_wd, "scenario3_f1_7methods_melt.rda"))


f1_facet_plot <- ggplot(f1_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    labs(title = NULL, x = NULL, y = "F1")

ggsave(
    file = paste0(figs_wd, "scenario3_F1_7methods_facet_plot.png"),
    plot = f1_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


rm(f1_melt, f1_facet_plot)

# --------------------
### POU
# --------------------
pou_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_fov", qfov, "_scenario3.rda")
        load(res_file)
        res1 <- t(sapply(res_anno_measures_scenario3$acc_pou, function(dat) dat["POU", paste0(methods, "unasg")]))
        res1 <- cbind(res_anno_measures_scenario3$acc_pou_iCAESAR["POU", "iCAESARunasg_ref1"], res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
pou_melt$Method <- factor(pou_melt$Method, levels = c("iCAESAR", methods))
pou_melt$Sample <- factor(pou_melt$Sample, levels = sample_used)

save(pou_melt, file = paste0(figs_wd, "scenario3_pou_7methods_melt.rda"))


pou_facet_plot <- ggplot(pou_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    labs(title = NULL, x = NULL, y = "POU")

ggsave(
    file = paste0(figs_wd, "scenario3_POU_7methods_facet_plot.png"),
    plot = pou_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


rm(pou_melt, pou_facet_plot)


# ==============================================================================
## control pou 0.05
# ==============================================================================

# --------------------
### ACC
# --------------------
acc_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_control_pou_0.05_fov", qfov, "_scenario3.rda")
        load(res_file)
        res1 <- t(sapply(res_anno_measures_scenario3$acc_pou, function(dat) dat["ACC", methods]))
        res1 <- cbind(res_anno_measures_scenario3$acc_pou_iCAESAR["ACC", "iCAESAR"], res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
acc_melt$Method <- factor(acc_melt$Method, levels = c("iCAESAR", methods))
acc_melt$Sample <- factor(acc_melt$Sample, levels = sample_used)

save(acc_melt, file = paste0(figs_wd, "scenario3_acc_control_pou_0.05_7methods_melt.rda"))



acc_facet_plot <- ggplot(acc_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    file = paste0(figs_wd, "scenario3_ACC_control_pou_0.05_7methods_facet_plot.png"),
    plot = acc_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


rm(acc_melt, acc_facet_plot)

# --------------------
### F1
# --------------------
f1_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_control_pou_0.05_fov", qfov, "_scenario3.rda")
        load(res_file)

        n_cts <- res_anno_measures_scenario3$n

        res1 <- t(sapply(res_anno_measures_scenario3$f1score, function(dat) {
            sapply(dat[methods], function(F1res) {
                nn <- n_cts[colnames(F1res)]
                nn <- nn / sum(nn)
                F1res["F1", ] %*% nn
            })
        }))

        F1res <- res_anno_measures_scenario3$f1score_iCAESAR[["iCAESAR"]]
        nn <- n_cts[colnames(F1res)]
        nn <- nn / sum(nn)
        F1_iCAESAR <- as.numeric(F1res["F1", ] %*% nn)
        
        res1 <- cbind(F1_iCAESAR, res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
f1_melt$Method <- factor(f1_melt$Method, levels = c("iCAESAR", methods))
f1_melt$Sample <- factor(f1_melt$Sample, levels = sample_used)

save(f1_melt, file = paste0(figs_wd, "scenario3_f1_control_pou_0.05_7methods_melt.rda"))


f1_facet_plot <- ggplot(f1_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    labs(title = NULL, x = NULL, y = "F1")

ggsave(
    file = paste0(figs_wd, "scenario3_F1_control_pou_0.05_7methods_facet_plot.png"),
    plot = f1_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


rm(f1_melt, f1_facet_plot)




# --------------------
### POU
# --------------------
pou_melt <- Reduce(rbind, lapply(2:4, function(i, methods) {
    # each query slice
    ana_wd <- paste0(work_wd, "processed_data_CosMx", i, "/")
    load(paste0(work_wd, "processed_data_CosMx", i, "/fovs.rda"))
    measures <- lapply(fovs, function(qfov) {
        res_file <- paste0(ana_wd, "fov", qfov, "/res_anno_measures_control_pou_0.05_fov", qfov, "_scenario3.rda")
        load(res_file)
        res1 <- t(sapply(res_anno_measures_scenario3$acc_pou, function(dat) dat["POU", methods]))
        res1 <- cbind(res_anno_measures_scenario3$acc_pou_iCAESAR["POU", "iCAESAR"], res1)
        res1
    })
    names(measures) <- fovs
    measures_array <- abind(measures, along = 3)
    dimnames(measures_array) <- list(
        rfov = 1:30,
        Method = c("iCAESAR", methods),
        qfov = fovs
    )
    res <- reshape2::melt(measures_array)
    cbind(res, Sample = sample_names[i])
}, methods = methods))
pou_melt$Method <- factor(pou_melt$Method, levels = c("iCAESAR", methods))
pou_melt$Sample <- factor(pou_melt$Sample, levels = sample_used)

save(pou_melt, file = paste0(figs_wd, "scenario3_pou_control_pou_0.05_7methods_melt.rda"))


pou_facet_plot <- ggplot(pou_melt, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
        position = position_dodge(width = 0.9), alpha = 0.5
    ) +
    geom_boxplot(
        outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.2
    ) +
    facet_wrap(~Sample, nrow = 1) +
    scale_fill_manual(labels = custom_labels, values = custom_colors) + # Apply custom colors
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
    labs(title = NULL, x = NULL, y = "POU")

ggsave(
    file = paste0(figs_wd, "scenario3_POU_control_pou_0.05_7methods_facet_plot.png"),
    plot = pou_facet_plot, width = 12.5, height = 2.5, units = "in", dpi = 200
)


rm(pou_melt, pou_facet_plot)





















