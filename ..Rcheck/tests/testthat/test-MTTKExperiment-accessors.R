test_that("replacement accessors update an MTTKExperiment", {
    counts <- matrix(
        c(1L, 2L, 3L, 4L),
        nrow = 2,
        dimnames = list(c("gene_1", "gene_2"), c("sample_1", "sample_2"))
    )

    x <- MTTKExperiment(assays = list(rna_counts = counts))

    genomeData(x) <- data.frame(
        genome_name = "Genome A",
        row.names = "genome_1"
    )
    links(x) <- list(
        gene_to_genome = data.frame(
            gene_id = c("gene_1", "gene_2"),
            genome_id = c("genome_1", "genome_1")
        )
    )
    activeHierarchies(x) <- "biological"
    dnaCounts(x) <- counts + 10L

    expect_identical(rownames(genomeData(x)), "genome_1")
    expect_identical(names(links(x)), "gene_to_genome")
    expect_identical(activeHierarchies(x), "biological")
    expect_identical(dim(dnaCounts(x)), c(2L, 2L))
    expect_true(methods::validObject(x))
})

test_that("replacement accessors can clear MTTK-specific metadata", {
    x <- make_test_mttk()

    genomeData(x) <- NULL
    links(x) <- NULL
    activeHierarchies(x) <- NULL

    expect_s4_class(genomeData(x), "DataFrame")
    expect_equal(nrow(genomeData(x)), 0L)
    expect_s4_class(links(x), "SimpleList")
    expect_length(links(x), 0L)
    expect_identical(activeHierarchies(x), character())
    expect_true(methods::validObject(x))
})

test_that("assay convenience accessors follow the named assay convention", {
    x <- make_test_mttk()

    expect_identical(
        names(SummarizedExperiment::assays(x, withDimnames = FALSE)),
        "rna_counts"
    )
    expect_null(dnaCounts(x))
    expect_null(rownames(rnaCounts(x, withDimnames = FALSE)))
    expect_null(colnames(rnaCounts(x, withDimnames = FALSE)))

    dnaCounts(x) <- rnaCounts(x) + 100L

    expect_identical(
        names(SummarizedExperiment::assays(x, withDimnames = FALSE)),
        c("rna_counts", "dna_counts")
    )
    expect_identical(dim(dnaCounts(x)), c(2L, 3L))
})
