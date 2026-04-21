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
        names(geneAssays(x, withDimnames = FALSE)),
        "rna_gene_counts"
    )
    expect_identical(
        names(genomeAssays(MTTKExample, withDimnames = FALSE)),
        "dna_genome_counts"
    )
    expect_identical(
        names(links(MTTKExample)),
        c(
            "gene_to_genome",
            "gene_to_ko",
            "ko_to_module",
            "ko_to_pathway",
            "module_to_pathway"
        )
    )
    expect_identical(rownames(genomeData(MTTKExample)), c("genome_1", "genome_2", "genome_3"))
    expect_identical(dim(dnaGenomeCounts(MTTKExample)), c(3L, 4L))
    expect_identical(S4Vectors::metadata(MTTKExample)$study_name,
        "Synthetic genome-resolved metatranscriptomics example"
    )
})

test_that("the packaged showcase object is available and reproducible", {
    x <- makeShowcaseMTTKExperiment()

    data("MTTKShowcase", envir = environment())

    expect_s4_class(x, "MTTKExperiment")
    expect_s4_class(MTTKShowcase, "MTTKExperiment")
    expect_identical(dim(x), c(16L, 6L))
    expect_identical(dim(MTTKShowcase), c(16L, 6L))
    expect_identical(names(links(x)), names(links(MTTKShowcase)))
    expect_identical(activeHierarchies(x), activeHierarchies(MTTKShowcase))
    expect_identical(
        names(geneAssays(x, withDimnames = FALSE)),
        "rna_gene_counts"
    )
    expect_identical(
        names(genomeAssays(MTTKShowcase, withDimnames = FALSE)),
        c("rna_genome_counts", "dna_genome_counts")
    )
    expect_identical(rownames(genomeData(MTTKShowcase)), c("genome_1", "genome_2", "genome_3", "genome_4"))
    expect_identical(dim(rnaGenomeCounts(MTTKShowcase)), c(4L, 6L))
    expect_identical(dim(dnaGenomeCounts(MTTKShowcase)), c(4L, 6L))
    domain_table <- table(as.character(genomeData(MTTKShowcase)$domain))
    expect_identical(names(domain_table), c("Archaea", "Bacteria"))
    expect_identical(as.integer(domain_table), c(2L, 2L))
    expect_identical(
        S4Vectors::metadata(MTTKShowcase)$study_name,
        "Synthetic coastal oxygen-pulse metatranscriptomics showcase"
    )
})
