
# 1. preprocess
source("MsHCC_analysis_1.R")

# 2. annotate and obtain fig. 3.d
for (i in 1:4) {
    commandArgs <- function(...) i
    source("MsHCC_analysis_2.R")
}

# 3. summary acc for fig. 3.e bottom panel
source("MsHCC_analysis_3_measures.R")

# 4. remove unwanted variation
source("MsHCC_analysis_4_ruv.R")

# 5. find integrate signature genes and obtain fig. 3.g
source("MsHCC_analysis_5_intsg.R")

# 6. obtain figure 3.f
source("MsHCC_analysis_6_plothcc.R")

# 7. obtain figure 3.h and figure 3.i
source("MsHCC_analysis_7_pathway.R")


















################################################################################
# combine acc to plot fig. 3.e
################################################################################
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)

load("hm_res_acc.rda")
load("ms_res_acc.rda")

myplot_embedding2 <- function(ACC, ylab) {
    data_df <- as.data.frame(ACC)
    data_df$Slice <- rownames(data_df)

    data_melted <- melt(data_df, id.vars = "Slice", variable.name = "Method", value.name = "Result")

    custom_colors <- c(
        "HCC1" = "#699ECA", "HCC2" = "#FF8C00",
        "HCC3" = "#F898CB", "HCC4" = "#4DAF4A"
    )
    custom_shapes <- c("HCC1" = 15, "HCC2" = 16, "HCC3" = 17, "HCC4" = 18)
    custom_labels <- c("CAESAR" = "CAESAR", "CELLID" = "Cell-ID")

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
            axis.title.y = element_text(face = "bold", size = 15), # Y-axis title
            legend.position = "bottom",
            legend.text = element_text(face = "bold", size = 12),
            legend.title = element_text(face = "bold", size = 12),
            axis.text.y = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 12)
        ) + # Position legend at the bottom
        labs(title = NULL, x = NULL, y = ylab)

    return(plot1)
}

plot1 <- myplot_embedding2(t(ms_res_acc), "ACC")
plot2 <- myplot_embedding2(t(hm_res_acc), "ACC")

library(patchwork)
combined_plot <- plot2 / plot1 + plot_layout(guides = "collect") & theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.height = unit(1, "lines"),
    legend.box.margin = margin(0, 0, 0, -15, unit = "mm")
) & guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(
    file = "acc.png",
    plot = combined_plot, width = 2.7, height = 5, units = "in", dpi = 200
)







