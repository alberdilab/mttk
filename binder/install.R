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
  "rmarkdown",
  "remotes"
), repos = "https://cloud.r-project.org")

repo_dir <- Sys.getenv("REPO_DIR", unset = ".")
remotes::install_local(repo_dir, dependencies = FALSE, upgrade = "never")
