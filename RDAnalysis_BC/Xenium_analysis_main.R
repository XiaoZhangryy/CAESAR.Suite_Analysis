
# 1. preprocess
source("Xenium_analysis_10_preprocess.R")

# 2. Fig 2.a 2.b
source("Xenium_analysis_11_Fig_2a_2b.R")

# 3. annotate and obtain fig. 3.c and fig. 3.e upper panel
for (i in 1:4) {
    commandArgs <- function(...) i
    source("Xenium_analysis_20_co-embedding_RCTD_annotation.R")
    
    commandArgs <- function(...) i
    source("Xenium_analysis_21_annotation_measures.R")
}

# 3. Fig 2.f
source("Xenium_analysis_22_Fig_2f_Supp_Data_1.R")

# 4. Fig 2.e
source("Xenium_analysis_23_Fig_2e.R")

# 5. RUV and Fig 2.c, Fig 2.d
source("Xenium_analysis_30_RUV.R")
source("Xenium_analysis_31_Fig2c.R")
source("Xenium_analysis_32_Fig_2d.R")

# 6. enrichment  and Fig 2.e, Fig 2.h
source("Xenium_analysis_40_run_enrichment.R")
source("Xenium_analysis_41_Fig_2e_2h.R")


