install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "BiocParallel",
  "SingleCellExperiment",
  "S4Vectors",
  "SummarizedExperiment",
  "TreeSummarizedExperiment",
  "BiocStyle",
  "KEGGREST"
), update = FALSE, ask = FALSE)

install.packages(c(
  "ape",
  "ggplot2",
  "glmmTMB",
  "nlme",
  "knitr",
  "rmarkdown"
), repos = "https://cloud.r-project.org")

