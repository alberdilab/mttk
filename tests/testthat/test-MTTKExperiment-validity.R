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

test_that("constructor requires authoritative genome ids in rowData", {
    components <- make_test_components()
    row_data <- components$rowData[, "gene_id", drop = FALSE]

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = row_data,
            colData = components$colData
        ),
        "rowData\\(x\\) must contain a 'genome_id' column"
    )
})

test_that("validity catches mismatched mirrored gene and sample identifiers", {
    components <- make_test_components()

    bad_row_data <- components$rowData
    bad_row_data$gene_id[2] <- "gene_x"

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = bad_row_data,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = components$links,
            activeHierarchies = components$activeHierarchies
        ),
        "rowData\\(x\\)\\$gene_id.*rownames\\(x\\)"
    )

    bad_col_data <- components$colData
    bad_col_data$sample_id[3] <- "sample_x"

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = bad_col_data,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = components$links,
            activeHierarchies = components$activeHierarchies
        ),
        "colData\\(x\\)\\$sample_id.*colnames\\(x\\)"
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

test_that("validity catches mismatches between rowData genome ids and gene_to_genome links", {
    components <- make_test_components()
    bad_links <- components$links
    bad_links$gene_to_genome$genome_id[2] <- "genome_1"

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = bad_links,
            activeHierarchies = components$activeHierarchies
        ),
        "must match rowData\\(x\\)\\$genome_id exactly"
    )
})

test_that("validity checks non-core link sources against canonical gene ids", {
    components <- make_test_components()
    bad_links <- components$links
    bad_links$gene_to_ko$gene_id[2] <- "gene_x"

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = bad_links,
            activeHierarchies = components$activeHierarchies
        ),
        "Link table 'gene_to_ko' contains gene_id values that are not defined"
    )
})

test_that("validity checks chained link namespaces against upstream targets", {
    components <- make_test_components()
    chained_links <- c(
        components$links,
        list(
            ko_to_module = data.frame(
                ko_id = c("K00001", "K99999"),
                module_id = c("M001", "M002")
            )
        )
    )

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = chained_links,
            activeHierarchies = components$activeHierarchies
        ),
        "Link table 'ko_to_module' contains ko_id values that are not defined"
    )
})

test_that("validity rejects unresolved source namespaces in link tables", {
    components <- make_test_components()
    bad_links <- c(
        components$links,
        list(
            orphan_to_group = data.frame(
                orphan_id = "orphan_a",
                group_id = "group_1"
            )
        )
    )

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = bad_links,
            activeHierarchies = components$activeHierarchies
        ),
        "unresolved source namespace 'orphan_id'"
    )
})

test_that("validity requires identifier column names in link tables", {
    components <- make_test_components()
    bad_links <- components$links
    bad_links$gene_to_ko <- data.frame(
        gene = c("gene_1", "gene_2"),
        ko_id = c("K00001", "K00002")
    )

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = components$genomeData,
            genomeAssays = components$genomeAssays,
            links = bad_links,
            activeHierarchies = components$activeHierarchies
        ),
        "must be identifier columns ending in '_id'"
    )
})

test_that("validity catches mismatched mirrored genome identifiers", {
    components <- make_test_components()
    genome_data <- S4Vectors::DataFrame(
        genome_id = c("genome_1", "genome_x"),
        genome_name = c("Genome A", "Genome B"),
        row.names = c("genome_1", "genome_2")
    )

    expect_error(
        MTTKExperiment(
            assays = components$assays,
            rowData = components$rowData,
            colData = components$colData,
            genomeData = genome_data,
            genomeAssays = components$genomeAssays,
            links = components$links,
            activeHierarchies = components$activeHierarchies
        ),
        "genomeData\\(x\\)\\$genome_id.*rownames\\(genomeExperiment\\(x\\)\\)"
    )
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
    row_data <- S4Vectors::DataFrame(
        genome_id = "genome_1",
        row.names = "gene_1"
    )
    col_data <- S4Vectors::DataFrame(
        row.names = "sample_1"
    )

    expect_error(
        MTTKExperiment(
            assays = list(rna_counts = counts),
            rowData = row_data,
            colData = col_data
        ),
        "explicit gene-level names"
    )

    x <- MTTKExperiment(
        assays = list(rna_gene_counts = counts),
        rowData = row_data,
        colData = col_data
    )
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
