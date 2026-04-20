test_that("MTTKExperiment constructor stores normalized metadata", {
    x <- make_test_mttk()

    expect_s4_class(x, "MTTKExperiment")
    expect_s4_class(geneExperiment(x), "TreeSummarizedExperiment")
    expect_s4_class(genomeExperiment(x), "SummarizedExperiment")
    expect_s4_class(genomeData(x), "DataFrame")
    expect_s4_class(links(x), "SimpleList")
    expect_identical(names(links(x)), c("gene_to_genome", "gene_to_ko"))
    expect_identical(activeHierarchies(x), c("biological", "functional"))
    expect_identical(dim(rnaGeneCounts(x)), c(2L, 3L))
    expect_identical(dim(dnaGenomeCounts(x)), c(2L, 3L))
})

test_that("constructor respects existing metadata(x)$mttk defaults", {
    metadata_input <- list(
        study_id = "test-study",
        mttk = list(
            links = methods::as(
                list(
                    gene_to_custom = S4Vectors::DataFrame(
                        gene_id = c("gene_1", "gene_2"),
                        custom_id = c("custom_a", "custom_b")
                    )
                ),
                "SimpleList"
            ),
            activeHierarchies = "custom_hierarchy"
        )
    )

    components <- make_test_components()
    x <- MTTKExperiment(
        assays = components$assays,
        rowData = components$rowData,
        colData = components$colData,
        metadata = metadata_input
    )

    expect_identical(S4Vectors::metadata(x)$study_id, "test-study")
    expect_null(genomeExperiment(x))
    expect_identical(nrow(genomeData(x)), 0L)
    expect_identical(names(links(x)), "gene_to_custom")
    expect_identical(activeHierarchies(x), "custom_hierarchy")
})

test_that("dedicated constructor arguments override metadata(x)$mttk", {
    metadata_input <- list(
        mttk = list(
            links = methods::as(
                list(
                    gene_to_custom = S4Vectors::DataFrame(
                        gene_id = c("gene_1", "gene_2"),
                        custom_id = c("custom_a", "custom_b")
                    )
                ),
                "SimpleList"
            ),
            activeHierarchies = "custom_hierarchy"
        )
    )

    x <- make_test_mttk(metadata = metadata_input)

    expect_identical(rownames(genomeData(x)), c("genome_1", "genome_2"))
    expect_identical(names(links(x)), c("gene_to_genome", "gene_to_ko"))
    expect_identical(activeHierarchies(x), c("biological", "functional"))
})

test_that("constructor rejects ambiguous assay names", {
    components <- make_test_components()
    bad_assays <- c(
        components$assays,
        list(dna_counts = components$assays$rna_gene_counts)
    )

    expect_error(
        MTTKExperiment(
            assays = bad_assays,
            rowData = components$rowData,
            colData = components$colData
        ),
        "Gene-level assays must use explicit gene-level names"
    )
})
