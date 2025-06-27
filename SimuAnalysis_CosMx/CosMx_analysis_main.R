
# 1. preprocess
for (i in 1:4) {
    commandArgs <- function(...) i
    source("CosMx_analysis_10_preprocess.R")
}

# 2. annotate and calculate measures for scenario 1
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_21_scenario1_annotation.R")
    
    commandArgs <- function(...) i
    source("CosMx_analysis_22_scenario1_annotation_measures.R")
}

# 3. summary measures and plot for scenario 1, obtain figure 1.f
source("CosMx_analysis_23_scenario1_Fig_1f.R")

# 4. pathway enrichment scoring
# markerList.rda used here as pathways for each cell type
for (i in 2:4) {
    commandArgs <- function(...) i
    source("CosMx_analysis_30_enrichment.R")
}

# 5. summary measures and plot for pathway enrichment scoring, obtain figure 1.d
source("CosMx_analysis_31_Fig_1g.R")

# 6. annotate and calculate measures for scenario 2
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_41_scenario2_annotation.R")
    
    commandArgs <- function(...) i
    source("CosMx_analysis_42_scenario2_annotation_measures.R")
}

# 7. summary measures and plot for scenario 2
source("CosMx_analysis_43_scenario2_plot_annotation_measures.R")

# 8. annotate and calculate measures for scenario 3
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_51_scenario3_annotation.R")
    
    commandArgs <- function(...) i
    source("CosMx_analysis_52_scenario3_annotation_measures.R")
}

# 9. summary measures and plot for scenario 3
source("CosMx_analysis_53_scenario3_plot_annotation_measures.R")






