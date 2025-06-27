rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/VisiumRDA/MsHCC/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

i <- as.integer(commandArgs(TRUE))

ana_wd <- paste0(work_wd, "HCC", i, "/")
setwd(ana_wd)

load("seu.rda")

################################################################################
# 1. CAESAR enrichment
################################################################################
load(paste0(work_wd, "pathwayList.rda"))

pathwaylist <- Reduce(c, pathwayList)

df_rgTest <- CAESAR.enrich.pathway(
    seu, pathway.list = pathwaylist, reduction = "caesar"
)
rownames(df_rgTest) <- names(pathwaylist)
save(df_rgTest, file = "df_rgTest.rda")

pathway_scores <- CAESAR.enrich.score(seu, pathwaylist)
save(pathway_scores, file = "pathway_scores.rda")
