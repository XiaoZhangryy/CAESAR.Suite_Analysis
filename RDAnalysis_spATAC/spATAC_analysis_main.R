
# 1. preprocess
source("spATAC_analysis_10_preprocess1.R")
source("spATAC_analysis_20_preprocess2.R")
source("spATAC_analysis_30_preprocess3.R")

# 2. annotate and obtain fig. 5.b
source("spATAC_analysis_40_coembedding_annotation_enrichscore_Fig_5b.R")

# 3. summary measures and obtain fig. 5.c
source("spATAC_analysis_41_annotation_measures.R")
source("spATAC_analysis_42_Fig_5c.R")

# 4. COUMAP and obtain fig. 5.a
source("spATAC_analysis_43_CoUMAP_Fig_5a.R")

# 5. enrichment and obtain fig. 5.d, 5.e
source("spATAC_analysis_50_enrichment_Fig_5d_5e.R")




