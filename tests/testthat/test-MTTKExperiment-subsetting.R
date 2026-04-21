test_that("row subsetting prunes genome-level rows and dependent links", {
    x <- makeExampleMTTKExperiment()

    y <- x[1:2, ]

    expect_s4_class(y, "MTTKExperiment")
    expect_identical(rownames(y), c("gene_1", "gene_2"))
    expect_identical(rownames(genomeData(y)), "genome_1")
    expect_identical(dim(dnaGenomeCounts(y)), c(1L, 4L))
    expect_identical(
        as.data.frame(links(y)$gene_to_genome),
        data.frame(
            gene_id = c("gene_1", "gene_2"),
            genome_id = c("genome_1", "genome_1"),
            row.names = c("gene_1", "gene_2")
        )
    )
    expect_identical(links(y)$gene_to_ko$gene_id, c("gene_1", "gene_2"))
    expect_identical(
        links(y)$ko_to_module$ko_id,
        c("K03043", "K03043", "K02111", "K02111")
    )
    expect_identical(
        links(y)$ko_to_pathway$ko_id,
        c("K03043", "K03043", "K02111", "K02111")
    )
    expect_identical(
        links(y)$module_to_pathway$module_id,
        c("M00001", "M00002", "M00009", "M00115")
    )
    expect_true(methods::validObject(y))
})

test_that("row subsetting also prunes the stored genome tree", {
    x <- makeShowcaseMTTKExperiment()

    y <- x[1:8, ]

    expect_s3_class(genomeTree(y), "phylo")
    expect_identical(sort(genomeTree(y)$tip.label), c("genome_1", "genome_2"))
    expect_identical(rownames(genomeData(y)), c("genome_1", "genome_2"))
    expect_true(methods::validObject(y))
})

test_that("row subsetting works without a stored gene_to_genome link table", {
    x <- makeExampleMTTKExperiment()
    kept_links <- links(x)
    kept_links$gene_to_genome <- NULL
    links(x) <- kept_links

    y <- x[3:4, ]

    expect_identical(rownames(y), c("gene_3", "gene_4"))
    expect_identical(rownames(genomeData(y)), "genome_2")
    expect_false("gene_to_genome" %in% names(links(y)))
    expect_identical(links(y)$gene_to_ko$gene_id, c("gene_3", "gene_4"))
    expect_true(methods::validObject(y))
})

test_that("column subsetting narrows both gene and genome assays", {
    x <- makeExampleMTTKExperiment()

    y <- x[, 1:2]

    expect_identical(dim(y), c(6L, 2L))
    expect_identical(dim(dnaGenomeCounts(y)), c(3L, 2L))
    expect_identical(rownames(genomeData(y)), c("genome_1", "genome_2", "genome_3"))
    expect_identical(names(links(y)), names(links(x)))
    expect_identical(links(y)$ko_to_pathway$ko_id, links(x)$ko_to_pathway$ko_id)
    expect_identical(links(y)$module_to_pathway$module_id, links(x)$module_to_pathway$module_id)
    expect_true(methods::validObject(y))
})

test_that("row subsetting to zero genes yields a valid empty nested object", {
    x <- makeExampleMTTKExperiment()

    y <- x[integer(), ]

    expect_identical(dim(y), c(0L, 4L))
    expect_identical(nrow(genomeData(y)), 0L)
    expect_identical(dim(dnaGenomeCounts(y)), c(0L, 4L))
    expect_identical(nrow(links(y)$gene_to_ko), 0L)
    expect_identical(nrow(links(y)$ko_to_module), 0L)
    expect_identical(nrow(links(y)$ko_to_pathway), 0L)
    expect_identical(nrow(links(y)$module_to_pathway), 0L)
    expect_true(methods::validObject(y))
})
