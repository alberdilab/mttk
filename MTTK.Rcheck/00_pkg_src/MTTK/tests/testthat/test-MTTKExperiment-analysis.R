test_that("aggregateToGenome sums assays and carries genome metadata", {
    x <- makeExampleMTTKExperiment()

    aggregated <- aggregateToGenome(x)

    expect_s4_class(aggregated, "SummarizedExperiment")
    expect_identical(dim(aggregated), c(3L, 4L))
    expect_identical(rownames(aggregated), c("genome_1", "genome_2", "genome_3"))
    expect_identical(
        names(SummarizedExperiment::assays(aggregated, withDimnames = FALSE)),
        "rna_genome_counts"
    )
    expect_identical(
        SummarizedExperiment::assay(aggregated, "rna_genome_counts"),
        rbind(
            genome_1 = rnaGeneCounts(x)["gene_1", ] + rnaGeneCounts(x)["gene_2", ],
            genome_2 = rnaGeneCounts(x)["gene_3", ] + rnaGeneCounts(x)["gene_4", ],
            genome_3 = rnaGeneCounts(x)["gene_5", ] + rnaGeneCounts(x)["gene_6", ]
        )
    )
    expect_true(all(c("genome_name", "clade", "taxonomy") %in%
        names(SummarizedExperiment::rowData(aggregated))))
    expect_identical(
        S4Vectors::metadata(aggregated)$mttk_aggregation$path,
        "gene_to_genome"
    )
})

test_that("aggregateByLink follows multi-step functional mappings", {
    x <- makeExampleMTTKExperiment()

    aggregated <- aggregateByLink(
        x,
        path = c("gene_to_ko", "ko_to_module"),
        assays = "rna_gene_counts"
    )

    expect_identical(rownames(aggregated), c("M001", "M002", "M003", "M004"))
    expect_identical(
        SummarizedExperiment::rowData(aggregated)$module_id,
        c("M001", "M002", "M003", "M004")
    )
    expect_identical(
        SummarizedExperiment::assay(aggregated, "rna_gene_counts"),
        rbind(
            M001 = rnaGeneCounts(x)["gene_1", ] + rnaGeneCounts(x)["gene_2", ],
            M002 = rnaGeneCounts(x)["gene_3", ],
            M003 = rnaGeneCounts(x)["gene_4", ],
            M004 = rnaGeneCounts(x)["gene_5", ] + rnaGeneCounts(x)["gene_6", ]
        )
    )
})

test_that("aggregateByLink can compute group means", {
    x <- makeExampleMTTKExperiment()

    aggregated <- aggregateByLink(
        x,
        path = "gene_to_genome",
        assays = "rna_gene_counts",
        fun = "mean"
    )

    expected <- rbind(
        genome_1 = (rnaGeneCounts(x)["gene_1", ] + rnaGeneCounts(x)["gene_2", ]) / 2,
        genome_2 = (rnaGeneCounts(x)["gene_3", ] + rnaGeneCounts(x)["gene_4", ]) / 2,
        genome_3 = (rnaGeneCounts(x)["gene_5", ] + rnaGeneCounts(x)["gene_6", ]) / 2
    )

    expect_equal(SummarizedExperiment::assay(aggregated, "rna_gene_counts"), expected)
})

test_that("summarizeActivity returns gene- and genome-level summaries", {
    x <- makeExampleMTTKExperiment()
    gene_to_genome <- SummarizedExperiment::rowData(x)$genome_id

    gene_activity <- summarizeActivity(x, by = "gene")
    genome_activity <- summarizeActivity(x, by = "genome")

    expect_s4_class(gene_activity, "SummarizedExperiment")
    expect_s4_class(genome_activity, "SummarizedExperiment")
    expect_identical(dim(gene_activity), c(6L, 4L))
    expect_identical(dim(genome_activity), c(3L, 4L))
    expect_identical(
        names(SummarizedExperiment::assays(gene_activity, withDimnames = FALSE)),
        "log2_activity"
    )
    expect_identical(
        names(SummarizedExperiment::assays(genome_activity, withDimnames = FALSE)),
        "log2_activity"
    )

    expected_gene <- log2(
        (rnaGeneCounts(x) + 1) /
            (dnaGenomeCounts(x)[match(gene_to_genome, rownames(dnaGenomeCounts(x))), , drop = FALSE] + 1)
    )
    rownames(expected_gene) <- rownames(x)

    expected_genome <- log2(
        (SummarizedExperiment::assay(
            aggregateToGenome(x, assays = "rna_gene_counts"),
            "rna_genome_counts"
        ) + 1) /
            (dnaGenomeCounts(x) + 1)
    )

    expect_equal(
        SummarizedExperiment::assay(gene_activity, "log2_activity"),
        expected_gene
    )
    expect_equal(
        SummarizedExperiment::assay(genome_activity, "log2_activity"),
        expected_genome
    )
    expect_identical(
        S4Vectors::metadata(genome_activity)$mttk_activity$by,
        "genome"
    )
})

test_that("analysis helpers reject unknown assays and links", {
    x <- makeExampleMTTKExperiment()

    expect_error(
        aggregateByLink(x, path = "missing_link"),
        "Unknown link table"
    )
    expect_error(
        aggregateByLink(x, path = "gene_to_genome", assays = "missing_assay"),
        "Unknown assay name"
    )
    expect_error(
        summarizeActivity(x, numeratorAssay = "missing_assay"),
        "Unknown assay name"
    )
    expect_error(
        summarizeActivity(x, genomeAssay = "missing_assay"),
        "Unknown genome assay name"
    )
})
