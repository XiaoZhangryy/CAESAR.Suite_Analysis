
# 1. preprocess
source("Xenium_analysis_1.R")

# 2. annotate and enrichment score
for (i in 1:4) {
    commandArgs <- function(...) i
    source("Xenium_analysis_2.R")
}

# 3. obtain figure 2.a and figure 2.b
source("Xenium_analysis_3_plotsc.R")

# 4. remove unwanted variation
source("Xenium_analysis_4_ruv.R")

# 5. summary measures and obtain supplementary data 1 and figure 2.f
source("Xenium_analysis_5_measures.R")

# 6. find integrate signature genes and obtain figure 2.d
source("Xenium_analysis_6_intsg.R")

# 7. obtain figure 2.c and figure 2.e
source("Xenium_analysis_7_plotbc.R")

# 8. obtain figure 2.g and figure 2.h
source("Xenium_analysis_8_pathway.R")











