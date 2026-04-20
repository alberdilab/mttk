test_that("aggregateToKOGenome returns KO/genome counts with stable metadata", {
    x <- makeExampleMTTKExperiment()
    aggregated <- aggregateToKOGenome(x)

    expect_s4_class(aggregated, "SummarizedExperiment")
    expect_identical(
        rownames(aggregated),
        c(
            "K03043::genome_1",
            "K02111::genome_1",
            "K03043::genome_2",
            "K00368::genome_2",
            "K02111::genome_3",
            "K00368::genome_3"
        )
    )
    expect_identical(
        as.character(SummarizedExperiment::rowData(aggregated)$ko_id),
        c("K03043", "K02111", "K03043", "K00368", "K02111", "K00368")
    )
    expect_identical(
        as.character(SummarizedExperiment::rowData(aggregated)$genome_id),
        c("genome_1", "genome_1", "genome_2", "genome_2", "genome_3", "genome_3")
    )
    expect_true(all(c("n_genes_pair", "n_links_pair") %in% names(SummarizedExperiment::rowData(aggregated))))
    expect_identical(
        rownames(S4Vectors::metadata(aggregated)$mttk_ko_aggregation$koSummary),
        c("K03043", "K02111", "K00368")
    )
    expect_equal(
        as.numeric(SummarizedExperiment::assay(aggregated, "rna_gene_counts")[1L, ]),
        as.numeric(rnaGeneCounts(x)["gene_1", ])
    )
})

test_that("makeKOMixedModelData aligns KO/genome observations and offsets", {
    x <- makeExampleMTTKExperiment()
    model_data <- makeKOMixedModelData(
        x,
        variable = "condition",
        model = "library_plus_genome_abundance"
    )

    expect_s4_class(model_data, "DataFrame")
    expect_identical(nrow(model_data), 24L)
    expect_true(all(c(
        "ko_genome_id",
        "sample_id",
        "rna_count",
        "ko_id",
        "genome_id",
        "condition",
        "lib_size",
        "genome_abundance",
        "genome_abundance_offset"
    ) %in% names(model_data)))

    model_state <- S4Vectors::metadata(model_data)$mttk_ko_model
    expect_identical(model_state$specification, "library_plus_genome_abundance")
    expect_identical(model_state$sourceAssay, "rna_gene_counts")
    expect_identical(model_state$libSizeSource, "colSums(rna_gene_counts)")
    expect_identical(model_state$genomeAssay, "dna_genome_counts")
    expect_identical(model_state$variableType, "two_level_factor")
    expect_identical(model_state$referenceLevel, "control")
    expect_identical(model_state$contrastLevel, "treated")
    expect_identical(model_state$effectLabel, "condition: treated vs control")

    one_row <- as.data.frame(
        model_data[
            model_data$ko_genome_id == "K03043::genome_1" &
                model_data$sample_id == "sample_1",
            ,
            drop = FALSE
        ]
    )
    expect_identical(one_row$ko_id, "K03043")
    expect_identical(one_row$genome_id, "genome_1")
    expect_equal(one_row$rna_count, 120)
    expect_equal(one_row$lib_size, sum(rnaGeneCounts(x)[, "sample_1"]))
    expect_equal(one_row$genome_abundance, dnaGenomeCounts(x)["genome_1", "sample_1"])
    expect_equal(one_row$genome_abundance_offset, one_row$genome_abundance + 1)
})

test_that("fitKOMixedModel fits the library-only KO mixed model", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKOMixedModel(x, variable = "condition")

    expect_s4_class(fit, "MTTKFit")
    expect_identical(rownames(fit), c("K03043", "K02111", "K00368"))
    expect_true(all(c(
        "ko_id",
        "effect_label",
        "estimate",
        "std_error",
        "statistic",
        "p_value",
        "q_value",
        "status",
        "n_genes",
        "n_links",
        "n_genomes"
    ) %in% names(fit)))
    expect_identical(fitInfo(fit)$backend, "glmmTMB")
    expect_identical(fitInfo(fit)$specification, "library_only")
    expect_identical(fitInfo(fit)$responseAssay, "rna_gene_counts")
    expect_identical(fitInfo(fit)$variableType, "two_level_factor")
    expect_identical(fitInfo(fit)$referenceLevel, "control")
    expect_identical(fitInfo(fit)$contrastLevel, "treated")
    expect_identical(fitInfo(fit)$effectLabel, "condition: treated vs control")
    expect_true(all(as.character(fit$effect_label) == "condition: treated vs control"))
    expect_length(modelObjects(fit), 0L)
    expect_true(all(as.character(fit$status) == "ok"))

    estimates <- stats::setNames(as.numeric(fit$estimate), rownames(fit))
    expect_gt(estimates[["K03043"]], 0)
    expect_lt(estimates[["K00368"]], 0)
    expect_true(all(!is.na(fit$q_value)))
})

test_that("fitKOMixedModel can include genome abundance offsets and keep models", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKOMixedModel(
        x,
        variable = "condition",
        model = "library_plus_genome_abundance",
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(
        fitInfo(fit)$specification,
        "library_plus_genome_abundance"
    )
    expect_identical(fitInfo(fit)$genomeAssay, "dna_genome_counts")
    expect_identical(fitInfo(fit)$offsetPseudocount, 1)
    expect_identical(sort(names(modelObjects(fit))), sort(rownames(fit)))
    expect_true(all(vapply(
        modelObjects(fit),
        function(x) methods::is(x, "glmmTMB"),
        logical(1)
    )))
})

test_that("fitKOMixedModel supports numeric sample-level variables", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKOMixedModel(x, variable = "pH")

    expect_s4_class(fit, "MTTKFit")
    expect_true(all(as.character(fit$status) == "ok"))
    expect_true(all(as.character(fit$tested_term) == "pH"))
    expect_identical(fitInfo(fit)$variableType, "numeric")
    expect_true(all(is.na(fit$reference_level)))
    expect_true(all(is.na(fit$contrast_level)))
    expect_true(all(as.character(fit$effect_label) == "pH per unit increase"))
})

test_that("fitKOMixedModel validates modeling inputs", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()

    expect_error(
        fitKOMixedModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitKOMixedModel(
            x,
            variable = "condition",
            model = "library_plus_genome_abundance",
            genomeAssay = "missing_assay"
        ),
        "Unknown genome assay name"
    )

    SummarizedExperiment::colData(x)$phase <- factor(c("a", "b", "c", "d"))
    expect_error(
        fitKOMixedModel(x, variable = "phase"),
        "factor with exactly two levels"
    )
})
