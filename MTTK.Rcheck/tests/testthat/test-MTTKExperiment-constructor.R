test_that("MTTKExperiment constructor stores normalized metadata", {
    x <- make_test_mttk()

    expect_s4_class(x, "MTTKExperiment")
    expect_s4_class(genomeData(x), "DataFrame")
    expect_s4_class(links(x), "SimpleList")
    expect_identical(names(links(x)), c("gene_to_genome", "gene_to_ko"))
    expect_identical(activeHierarchies(x), c("biological", "functional"))
    expect_identical(dim(rnaCounts(x)), c(2L, 3L))
    expect_null(dnaCounts(x))
})

test_that("constructor respects existing metadata(x)$mttk defaults", {
    metadata_input <- list(
        study_id = "test-study",
        mttk = list(
            genomeData = S4Vectors::DataFrame(
                genome_name = "Genome X",
                row.names = "genome_x"
            ),
            links = methods::as(
                list(
                    gene_to_genome = S4Vectors::DataFrame(
                        gene_id = c("gene_1", "gene_2"),
                        genome_id = c("genome_x", "genome_x")
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
    expect_identical(rownames(genomeData(x)), "genome_x")
    expect_identical(names(links(x)), "gene_to_genome")
    expect_identical(activeHierarchies(x), "custom_hierarchy")
})

test_that("dedicated constructor arguments override metadata(x)$mttk", {
    metadata_input <- list(
        mttk = list(
            genomeData = S4Vectors::DataFrame(
                genome_name = "Genome X",
                row.names = "genome_x"
            ),
            links = methods::as(
                list(
                    gene_to_genome = S4Vectors::DataFrame(
                        gene_id = c("gene_1", "gene_2"),
                        genome_id = c("genome_x", "genome_x")
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
