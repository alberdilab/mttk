test_that("links must be named on construction", {
    counts <- matrix(1L, nrow = 1, ncol = 1)

    expect_error(
        MTTKExperiment(
            assays = list(rna_counts = counts),
            links = list(data.frame(gene_id = "gene_1", genome_id = "genome_1"))
        ),
        "named collection"
    )
})

test_that("replacement accessors reject invalid values", {
    x <- make_test_mttk()

    expect_error(
        genomeData(x) <- matrix(1:4, nrow = 2),
        "must be an S4Vectors::DataFrame, a data.frame, or NULL"
    )
    expect_error(
        links(x) <- list(gene_to_genome = matrix(1:4, nrow = 2)),
        "must be an S4Vectors::DataFrame or data.frame"
    )
    expect_error(
        activeHierarchies(x) <- c("biological", ""),
        "must contain non-missing, non-empty names"
    )
})

test_that("validity catches malformed metadata(x)$mttk values", {
    x <- make_test_mttk()
    metadata_list <- S4Vectors::metadata(x)

    metadata_list$mttk <- "not_a_list"
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "must be a list")

    metadata_list <- S4Vectors::metadata(x)
    bad_genome_data <- S4Vectors::DataFrame(genome_name = c("Genome A", "Genome B"))
    rownames(bad_genome_data) <- c("genome_1", "genome_1")
    metadata_list$mttk <- list(
        genomeData = bad_genome_data,
        links = S4Vectors::SimpleList(),
        activeHierarchies = character()
    )
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "Row names of 'genomeData' must be unique")

    metadata_list <- S4Vectors::metadata(x)
    metadata_list$mttk <- list(
        genomeData = S4Vectors::DataFrame(),
        links = S4Vectors::SimpleList(gene_to_genome = matrix(1:4, nrow = 2)),
        activeHierarchies = character()
    )
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "Each entry in 'links' must be an S4Vectors::DataFrame")

    metadata_list <- S4Vectors::metadata(x)
    metadata_list$mttk <- list(
        genomeData = S4Vectors::DataFrame(),
        links = S4Vectors::SimpleList(),
        activeHierarchies = c("biological", "biological")
    )
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "'activeHierarchies' must contain unique, non-empty names")
})
