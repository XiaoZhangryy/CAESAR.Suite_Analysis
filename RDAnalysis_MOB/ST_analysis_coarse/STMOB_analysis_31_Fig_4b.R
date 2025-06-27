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


plot_methods <- c(
  "CAESAR", "CelliD", "Seurat", "scmap", "SingleR", "scPred"
)
custom_labels <- custom_labels[plot_methods]
custom_colors <- custom_colors[plot_methods]


#===============================================================================
# 1.1 plot measures - without account for unassigned
#===============================================================================
load("res_acc_pou_f1_7methods.rda")

## --------------------
## ACC
## --------------------

measures <- t(res_acc_pou["ACC", methods])
dimnames(measures) <- list(
    Sample = "ST",
    Method = methods
)
acc_melt <- reshape2::melt(measures)
acc_melt$Method <- factor(acc_melt$Method, levels = methods)

save(acc_melt, file = "acc_7methods_melt.rda")



load("acc_7methods_melt.rda")

acc_melt <- acc_melt[acc_melt$Method %in% plot_methods, ]
acc_melt$Method <- factor(acc_melt$Method, levels = plot_methods)

acc_facet_plot <- ggplot(acc_melt, aes(x = Method, y = value, fill = Method)) +
  geom_bar(
    stat = "identity", width = 0.8, color = "black"
  ) +
  geom_text(aes(label = round(value, 3)),
          vjust = -0.5,
          size = 4,
          fontface = "bold") +
  scale_fill_manual(labels = custom_labels, values = custom_colors) +
  scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10, angle = 25, hjust = 1),
    strip.text.x = element_text(size = 10)
  ) +
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = "ACC")


ggsave(
    file = paste0(figs_wd, "ACC_7methods_facet_plot.png"),
    plot = acc_facet_plot, width = 6, height = 4, units = "in", dpi = 200
)

ggsave(
    file = paste0(figs_wd, "ACC_7methods_facet_plot_nolegend.png"),
    plot = acc_facet_plot + theme(legend.position = "none"), width = 5, height = 4, units = "in", dpi = 200
)



## --------------------
## F1
## --------------------

measures <- t(f1score[methods])
dimnames(measures) <- list(
    Sample = "ST",
    Method = methods
)
f1_melt <- reshape2::melt(measures)
f1_melt$Method <- factor(f1_melt$Method, levels = methods)

save(f1_melt, file = "f1_7methods_melt.rda")



load("f1_7methods_melt.rda")

f1_melt <- f1_melt[f1_melt$Method %in% plot_methods, ]
f1_melt$Method <- factor(f1_melt$Method, levels = plot_methods)

f1_facet_plot <- ggplot(f1_melt, aes(x = Method, y = value, fill = Method)) +
  geom_bar(
    stat = "identity", width = 0.8, color = "black"
  ) +
  geom_text(aes(label = round(value, 3)),
          vjust = -0.5,
          size = 4,
          fontface = "bold") +
  scale_fill_manual(labels = custom_labels, values = custom_colors) +
  scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10, angle = 25, hjust = 1),
    strip.text.x = element_text(size = 10)
  ) +
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = "F1")


ggsave(
    file = paste0(figs_wd, "F1_7methods_facet_plot.png"),
    plot = f1_facet_plot, width = 6, height = 4, units = "in", dpi = 200
)

ggsave(
    file = paste0(figs_wd, "F1_7methods_facet_plot_nolegend.png"),
    plot = f1_facet_plot + theme(legend.position = "none"), width = 5, height = 4, units = "in", dpi = 200
)




#===============================================================================
# 1.3 plot measures - control pou 0.05
#===============================================================================
load("res_acc_pou_f1_7methods_control_pou_0.05.rda")

## --------------------
## ACC
## --------------------

measures <- t(res_acc_pou["ACC", methods])
dimnames(measures) <- list(
    Sample = "ST",
    Method = methods
)
acc_melt <- reshape2::melt(measures)
acc_melt$Method <- factor(acc_melt$Method, levels = methods)

save(acc_melt, file = "acc_control_pou_0.05_7methods_melt.rda")



load("acc_control_pou_0.05_7methods_melt.rda")

acc_melt <- acc_melt[acc_melt$Method %in% plot_methods, ]
acc_melt$Method <- factor(acc_melt$Method, levels = plot_methods)

acc_facet_plot <- ggplot(acc_melt, aes(x = Method, y = value, fill = Method)) +
  geom_bar(
    stat = "identity", width = 0.8, color = "black"
  ) +
  geom_text(aes(label = round(value, 3)),
          vjust = -0.5,
          size = 4,
          fontface = "bold") +
  scale_fill_manual(labels = custom_labels, values = custom_colors) +
  scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10, angle = 25, hjust = 1),
    strip.text.x = element_text(size = 10)
  ) +
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = "ACC")


ggsave(
    file = paste0(figs_wd, "ACC_control_pou_0.05_7methods_facet_plot.png"),
    plot = acc_facet_plot, width = 6, height = 4, units = "in", dpi = 200
)

ggsave(
    file = paste0(figs_wd, "ACC_control_pou_0.05_7methods_facet_plot_nolegend.png"),
    plot = acc_facet_plot + theme(legend.position = "none"), width = 5, height = 4, units = "in", dpi = 200
)



## --------------------
## F1
## --------------------

measures <- t(f1score[methods])
dimnames(measures) <- list(
    Sample = "ST",
    Method = methods
)
f1_melt <- reshape2::melt(measures)
f1_melt$Method <- factor(f1_melt$Method, levels = methods)

save(f1_melt, file = "f1_control_pou_0.05_7methods_melt.rda")



load("f1_control_pou_0.05_7methods_melt.rda")

f1_melt <- f1_melt[f1_melt$Method %in% plot_methods, ]
f1_melt$Method <- factor(f1_melt$Method, levels = plot_methods)

f1_facet_plot <- ggplot(f1_melt, aes(x = Method, y = value, fill = Method)) +
  geom_bar(
    stat = "identity", width = 0.8, color = "black"
  ) +
  geom_text(aes(label = round(value, 3)),
          vjust = -0.5,
          size = 4,
          fontface = "bold") +
  scale_fill_manual(labels = custom_labels, values = custom_colors) +
  scale_x_discrete(labels = custom_labels, limits = names(custom_labels)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10, angle = 25, hjust = 1),
    strip.text.x = element_text(size = 10)
  ) +
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = "F1")


ggsave(
    file = paste0(figs_wd, "F1_control_pou_0.05_7methods_facet_plot.png"),
    plot = f1_facet_plot, width = 6, height = 4, units = "in", dpi = 200
)

ggsave(
    file = paste0(figs_wd, "F1_control_pou_0.05_7methods_facet_plot_nolegend.png"),
    plot = f1_facet_plot + theme(legend.position = "none"), width = 5, height = 4, units = "in", dpi = 200
)

