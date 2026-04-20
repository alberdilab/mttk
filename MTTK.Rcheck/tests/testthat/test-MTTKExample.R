test_that("the packaged example object is available and reproducible", {
    x <- makeExampleMTTKExperiment()

    data("MTTKExample", envir = environment())

    expect_s4_class(x, "MTTKExperiment")
    expect_s4_class(MTTKExample, "MTTKExperiment")
    expect_identical(dim(x), c(6L, 4L))
    expect_identical(dim(MTTKExample), c(6L, 4L))
    expect_identical(names(links(x)), names(links(MTTKExample)))
    expect_identical(activeHierarchies(x), activeHierarchies(MTTKExample))
    expect_identical(
        names(SummarizedExperiment::assays(x, withDimnames = FALSE)),
        c("rna_counts", "dna_counts")
    )
    expect_identical(
        names(links(MTTKExample)),
        c("gene_to_genome", "gene_to_ko", "ko_to_module", "module_to_pathway")
    )
    expect_identical(rownames(genomeData(MTTKExample)), c("genome_1", "genome_2", "genome_3"))
    expect_identical(S4Vectors::metadata(MTTKExample)$study_name,
        "Synthetic genome-resolved metatranscriptomics example"
    )
})
