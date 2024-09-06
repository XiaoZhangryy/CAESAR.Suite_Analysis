
# 1. preprocess
source("CosMx_analysis_1.R")

# 2. annotate and calculate measures for scenario 1
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_2.R")
}

# 3. summary measures and plot for scenario 1, obtain figure 1.c left penal
source("CosMx_analysis_3.R")

# 4. pathway enrichment scoring
for (i in 2:4) {
    commandArgs <- function(...) i
    source("CosMx_analysis_4.R")
}

# 5. summary measures and plot for pathway enrichment scoring, obtain figure 1.d
source("CosMx_analysis_5.R")

# 6. annotate and calculate measures for scenario 2
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_6.R")
}

# 7. summary measures and plot for scenario 2, obtain figure 1.c right penal
source("CosMx_analysis_7.R")

# 8. annotate and calculate measures for scenario 3
for (i in 1:77) {
    commandArgs <- function(...) i
    source("CosMx_analysis_8.R")
}

# 9. summary measures and plot for scenario 3
source("CosMx_analysis_9.R")






