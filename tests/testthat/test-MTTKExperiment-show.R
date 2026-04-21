test_that("show method summarizes the main MTTKExperiment components", {
    x <- make_test_mttk()

    output <- capture.output(show(x))

    expect_true(any(grepl("^MTTKExperiment with 2 features and 3 samples$", output)))
    expect_true(any(grepl("^geneAssays\\(1\\): rna_gene_counts$", output)))
    expect_true(any(grepl("^rowData names\\(2\\): gene_id, genome_id$", output)))
    expect_true(any(grepl("^colData names\\(2\\): sample_id, condition$", output)))
    expect_true(any(grepl("^genomeData rows\\(2\\): 2$", output)))
    expect_true(any(grepl("^genomeAssays\\(1\\): dna_genome_counts$", output)))
    expect_true(any(grepl("^genomeTree tips\\(0\\): 0$", output)))
    expect_true(any(grepl("^links\\(2\\): gene_to_genome, gene_to_ko$", output)))
    expect_true(any(grepl("^activeHierarchies\\(2\\): biological, functional$", output)))
    expect_true(any(grepl("^rowTreeNames\\(0\\): none$", output)))
    expect_true(any(grepl("^colTreeNames\\(0\\): none$", output)))
})
