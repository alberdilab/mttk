test_that("links must be named on construction", {
    counts <- matrix(1L, nrow = 1, ncol = 1)

    expect_error(
        MTTKExperiment(
            assays = list(rna_gene_counts = counts),
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
    metadata_list$mttk <- list(
        links = S4Vectors::SimpleList(gene_to_genome = matrix(1:4, nrow = 2)),
        activeHierarchies = character()
    )
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "Each entry in 'links' must be an S4Vectors::DataFrame")

    metadata_list <- S4Vectors::metadata(x)
    metadata_list$mttk <- list(
        links = S4Vectors::SimpleList(),
        activeHierarchies = c("biological", "biological")
    )
    S4Vectors::metadata(x) <- metadata_list
    expect_error(methods::validObject(x), "'activeHierarchies' must contain unique, non-empty names")
})

test_that("validity catches genome mappings that are missing from genomeExperiment", {
    components <- make_test_components()
    bad_genome_data <- components$genomeData["genome_1", , drop = FALSE]
    bad_genome_assays <- list(
        dna_genome_counts = components$genomeAssays$dna_genome_counts["genome_1", , drop = FALSE]
    )

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = bad_genome_data,
            genomeAssays = bad_genome_assays,
            links = components$links,
            activeHierarchies = components$activeHierarchies
        ),
        "All mapped genomes must be present in 'genomeExperiment\\(x\\)'"
    )
})

test_that("validity rejects ambiguous assay names", {
    counts <- matrix(1L, nrow = 1, ncol = 1, dimnames = list("gene_1", "sample_1"))

    expect_error(
        MTTKExperiment(assays = list(rna_counts = counts)),
        "explicit gene-level names"
    )

    x <- MTTKExperiment(assays = list(rna_gene_counts = counts))
    expect_error(
        genomeExperiment(x) <- SummarizedExperiment::SummarizedExperiment(
            assays = list(dna_counts = matrix(
                1L,
                nrow = 1,
                ncol = 1,
                dimnames = list("genome_1", "sample_1")
            ))
        ),
        "explicit genome-level names"
    )
})
