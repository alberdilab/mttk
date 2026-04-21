test_that("replacement accessors update an MTTKExperiment", {
    gene_counts <- matrix(
        c(1L, 2L, 3L, 4L),
        nrow = 2,
        dimnames = list(c("gene_1", "gene_2"), c("sample_1", "sample_2"))
    )
    dna_counts <- matrix(
        c(11L, 12L),
        nrow = 1,
        dimnames = list("genome_1", c("sample_1", "sample_2"))
    )
    row_data <- S4Vectors::DataFrame(
        genome_id = c("genome_1", "genome_1"),
        row.names = c("gene_1", "gene_2")
    )
    col_data <- S4Vectors::DataFrame(
        row.names = c("sample_1", "sample_2")
    )

    x <- MTTKExperiment(
        assays = list(rna_gene_counts = gene_counts),
        rowData = row_data,
        colData = col_data
    )

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
    dnaGenomeCounts(x) <- dna_counts

    expect_identical(rownames(genomeData(x)), "genome_1")
    expect_identical(names(links(x)), "gene_to_genome")
    expect_identical(activeHierarchies(x), "biological")
    expect_identical(dim(dnaGenomeCounts(x)), c(1L, 2L))
    expect_identical(
        names(geneAssays(x, withDimnames = FALSE)),
        "rna_gene_counts"
    )
    expect_true(methods::validObject(x))
})

test_that("replacement accessors can clear MTTK-specific metadata", {
    x <- make_test_mttk()

    genomeData(x) <- NULL
    links(x) <- NULL
    activeHierarchies(x) <- NULL

    expect_s4_class(genomeData(x), "DataFrame")
    expect_equal(nrow(genomeData(x)), 2L)
    expect_equal(ncol(genomeData(x)), 0L)
    expect_s4_class(links(x), "SimpleList")
    expect_length(links(x), 0L)
    expect_identical(activeHierarchies(x), character())
    expect_true(methods::validObject(x))
})

test_that("declared gene and genome accessors expose the two data layers", {
    x <- make_test_mttk()
    raw_counts <- rnaGeneCounts(x, withDimnames = FALSE)

    expect_identical(
        names(geneAssays(x, withDimnames = FALSE)),
        "rna_gene_counts"
    )
    expect_identical(
        names(genomeAssays(x, withDimnames = FALSE)),
        "dna_genome_counts"
    )
    expect_identical(dim(dnaGenomeCounts(x)), c(2L, 3L))
    expect_identical(dim(raw_counts), c(2L, 3L))
    expect_identical(unname(raw_counts), unname(rnaGeneCounts(x)))
    expect_s4_class(geneExperiment(x), "TreeSummarizedExperiment")
    expect_s4_class(genomeExperiment(x), "SummarizedExperiment")

    dnaGenomeCounts(x) <- dnaGenomeCounts(x) + 100L

    expect_identical(
        names(geneAssays(x, withDimnames = FALSE)),
        "rna_gene_counts"
    )
    expect_identical(
        names(genomeAssays(x, withDimnames = FALSE)),
        "dna_genome_counts"
    )
    expect_identical(dim(dnaGenomeCounts(x)), c(2L, 3L))
})

test_that("genomeTree can be set, retrieved, and create a genome layer when needed", {
    x <- MTTKExperiment(
        assays = make_test_components()$assays,
        rowData = make_test_components()$rowData,
        colData = make_test_components()$colData
    )
    tree <- ape::read.tree(text = "(genome_1:0.1,genome_2:0.1)root;")

    genomeTree(x) <- tree

    expect_s4_class(genomeExperiment(x), "SummarizedExperiment")
    expect_s3_class(genomeTree(x), "phylo")
    expect_identical(sort(genomeTree(x)$tip.label), c("genome_1", "genome_2"))
    expect_identical(rownames(genomeData(x)), c("genome_1", "genome_2"))
    expect_true(methods::validObject(x))
})

test_that("genome-level count accessors can create and update the genome experiment", {
    counts <- matrix(
        c(1L, 2L, 3L, 4L),
        nrow = 2,
        dimnames = list(c("gene_1", "gene_2"), c("sample_1", "sample_2"))
    )
    genome_counts <- matrix(
        c(20L, 25L),
        nrow = 1,
        dimnames = list("genome_1", c("sample_1", "sample_2"))
    )
    row_data <- S4Vectors::DataFrame(
        genome_id = c("genome_1", "genome_1"),
        row.names = c("gene_1", "gene_2")
    )
    col_data <- S4Vectors::DataFrame(
        row.names = c("sample_1", "sample_2")
    )

    x <- MTTKExperiment(
        assays = list(rna_gene_counts = counts),
        rowData = row_data,
        colData = col_data
    )
    dnaGenomeCounts(x) <- genome_counts
    rnaGenomeCounts(x) <- genome_counts + 100L

    expect_s4_class(genomeExperiment(x), "SummarizedExperiment")
    expect_identical(rownames(genomeData(x)), "genome_1")
    expect_identical(dnaGenomeCounts(x), genome_counts)
    expect_identical(rnaGenomeCounts(x), genome_counts + 100L)
})
