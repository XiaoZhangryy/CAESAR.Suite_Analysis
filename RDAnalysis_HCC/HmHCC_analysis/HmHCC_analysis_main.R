
# 1. preprocess
source("HmHCC_analysis_1.R")

# 2. annotate and obtain fig. 3.c and fig. 3.e upper panel
for (i in 1:4) {
    commandArgs <- function(...) i
    source("HmHCC_analysis_2.R")

    commandArgs <- function(...) i
    source("HmHCC_analysis_21_annotation_measures.R")
}

source("HmHCC_analysis_22_Fig_3e.R")





