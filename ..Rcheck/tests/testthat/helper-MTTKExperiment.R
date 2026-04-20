make_test_components <- function() {
    counts <- matrix(
        c(10L, 12L, 3L, 4L, 20L, 21L),
        nrow = 2,
        dimnames = list(c("gene_1", "gene_2"), c("sample_1", "sample_2", "sample_3"))
    )

    row_data <- S4Vectors::DataFrame(
        gene_id = rownames(counts),
        genome_id = c("genome_1", "genome_2"),
        row.names = rownames(counts)
    )

    col_data <- S4Vectors::DataFrame(
        sample_id = colnames(counts),
        condition = c("control", "control", "treated"),
        row.names = colnames(counts)
    )

    genome_data <- data.frame(
        genome_name = c("Genome A", "Genome B"),
        row.names = c("genome_1", "genome_2")
    )

    link_tables <- list(
        gene_to_genome = data.frame(
            gene_id = rownames(counts),
            genome_id = c("genome_1", "genome_2")
        ),
        gene_to_ko = data.frame(
            gene_id = rownames(counts),
            ko_id = c("K00001", "K00002")
        )
    )

    list(
        assays = list(rna_counts = counts),
        rowData = row_data,
        colData = col_data,
        genomeData = genome_data,
        links = link_tables,
        activeHierarchies = c("biological", "functional")
    )
}

make_test_mttk <- function(
    assays = NULL,
    genomeData = NULL,
    links = NULL,
    activeHierarchies = NULL,
    metadata = list()
) {
    components <- make_test_components()

    if (is.null(assays)) {
        assays <- components$assays
    }

    if (is.null(genomeData)) {
        genomeData <- components$genomeData
    }

    if (is.null(links)) {
        links <- components$links
    }

    if (is.null(activeHierarchies)) {
        activeHierarchies <- components$activeHierarchies
    }

    MTTKExperiment(
        assays = assays,
        rowData = components$rowData,
        colData = components$colData,
        genomeData = genomeData,
        links = links,
        activeHierarchies = activeHierarchies,
        metadata = metadata
    )
}
