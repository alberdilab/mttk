if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to build example data.", call. = FALSE)
}

pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

MTTKExample <- makeExampleMTTKExperiment()

if (!dir.exists("data")) {
    dir.create("data", recursive = TRUE)
}

save(MTTKExample, file = "data/MTTKExample.rda", compress = "xz")
