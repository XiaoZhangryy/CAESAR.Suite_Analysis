# CAESAR: a cross-technology and cross-resolution framework for spatial omics annotation

## Validation using CosMx data
Brief descriptions of simulated scripts (./SimuAnalysis_CosMx folder):

**CosMx_analysis_10_preprocess.R**: Preprocessing of CosMx data for simulation.

**CosMx_analysis_21_scenario1_annotation.R**, **CosMx_analysis_22_scenario1_annotation_measures.R** and **CosMx_analysis_23_scenario1_Fig_1f.R**: Simulation analysis for Scenario 1 with high-resolution reference and high-resolution target data.

**CosMx_analysis_30_enrichment.R** and **CosMx_analysis_31_Fig_1g.R**: Simulation analysis for pathway enrichment scoring with high-resolution reference and high-resolution target data.

**CosMx_analysis_41_scenario2_annotation.R**, **CosMx_analysis_42_scenario2_annotation_measures.R** and **CosMx_analysis_43_scenario2_Fig_1f.R**: Simulation analysis for Scenario 2 with high-resolution reference and low-resolution target data.

**CosMx_analysis_51_scenario3_annotation.R**, **CosMx_analysis_52_scenario3_annotation_measures.R** and **CosMx_analysis_53_scenario3_Fig_1f.R**: Simulation analysis for Scenario 3 with low-resolution reference and high-resolution target data.

**CosMx_analysis_main.R**: Process of running simulation.

## Real data analysis

**RDAnalysis_BC/Xenium_analysis_main.R**: Analysis for breast cancer Xenium data. 

**RDAnalysis_HCC/HmHCC_analysis/HmHCC_analysis_main.R**: Analysis for human hepatocellular carcinoma Visium data using human scRNA-seq reference. 

**RDAnalysis_HCC/MsHCC_analysis/MsHCC_analysis_main.R**: Analysis for human hepatocellular carcinoma Visium data using mouse scRNA-seq reference. 

**RDAnalysis_MOB/ST_analysis_coarse/STMOB_analysis_main.R**: Analysis for mouse olfactory bulb ST data using coarse-grained labelled scRNA-seq reference. 

**RDAnalysis_MOB/Pixel_analysis_fine/PixelMOB_analysis_main.R**: Analysis for mouse olfactory bulb Pixel-seq data using fine-grained labelled scRNA-seq reference. 

**RDAnalysis_spATAC/spATAC_analysis_main.R**: Analysis for mouse embryo E11 spatial ATAC-seq data. 

## Real data results 
The real data outputs (RD_results folder).

**Supplementray Data 1.xlsx**: Marker list for breast cancer Xenium data.

