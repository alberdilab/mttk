test_that("makeGeneModelData aligns gene/sample observations and offsets", {
    x <- makeExampleMTTKExperiment()
    model_data <- makeGeneModelData(
        x,
        variable = "condition",
        genomeOffset = TRUE
    )

    expect_s4_class(model_data, "DataFrame")
    expect_identical(nrow(model_data), 24L)
    expect_true(all(c(
        "gene_id",
        "sample_id",
        "rna_count",
        "genome_id",
        "gene_name",
        "condition",
        "lib_size",
        "genome_abundance",
        "genome_abundance_offset"
    ) %in% names(model_data)))

    model_state <- S4Vectors::metadata(model_data)$mttk_gene_model
    expect_identical(model_state$specification, "library_plus_genome_abundance")
    expect_identical(model_state$libraryOffset, TRUE)
    expect_identical(model_state$genomeOffset, TRUE)
    expect_identical(model_state$sourceAssay, "rna_gene_counts")
    expect_identical(model_state$libSizeSource, "colSums(rna_gene_counts)")
    expect_identical(model_state$genomeAssay, "dna_genome_counts")
    expect_identical(model_state$variableType, "two_level_factor")
    expect_identical(model_state$referenceLevel, "control")
    expect_identical(model_state$contrastLevel, "treated")
    expect_identical(model_state$effectLabel, "condition: treated vs control")

    one_row <- as.data.frame(
        model_data[
            model_data$gene_id == "gene_1" &
                model_data$sample_id == "sample_1",
            ,
            drop = FALSE
        ]
    )
    expect_identical(one_row$gene_id, "gene_1")
    expect_identical(one_row$genome_id, "genome_1")
    expect_identical(one_row$gene_name, "rpoB_like")
    expect_equal(one_row$rna_count, 120)
    expect_equal(one_row$lib_size, sum(rnaGeneCounts(x)[, "sample_1"]))
    expect_equal(one_row$genome_abundance, dnaGenomeCounts(x)["genome_1", "sample_1"])
    expect_equal(one_row$genome_abundance_offset, one_row$genome_abundance + 1)
})

test_that("fitGeneModel fits the library-only gene model", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitGeneModel(x, variable = "condition", genomeOffset = FALSE)

    expect_s4_class(fit, "MTTKFit")
    expect_identical(rownames(fit), paste0("gene_", seq_len(6)))
    expect_true(all(c(
        "gene_id",
        "genome_id",
        "gene_name",
        "effect_label",
        "estimate",
        "std_error",
        "statistic",
        "p_value",
        "q_value",
        "status"
    ) %in% names(fit)))
    expect_identical(fitInfo(fit)$backend, "glmmTMB")
    expect_identical(fitInfo(fit)$model, "gene_model")
    expect_identical(fitInfo(fit)$specification, "library_only")
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, FALSE)
    expect_identical(fitInfo(fit)$responseAssay, "rna_gene_counts")
    expect_identical(fitInfo(fit)$variableType, "two_level_factor")
    expect_identical(fitInfo(fit)$referenceLevel, "control")
    expect_identical(fitInfo(fit)$contrastLevel, "treated")
    expect_identical(fitInfo(fit)$effectLabel, "condition: treated vs control")
    expect_true(all(as.character(fit$effect_label) == "condition: treated vs control"))
    expect_true(all(as.character(fit$status) == "ok"))
    expect_true(all(!is.na(fit$q_value)))

    estimates <- stats::setNames(as.numeric(fit$estimate), rownames(fit))
    expect_gt(estimates[["gene_1"]], 0)
    expect_lt(estimates[["gene_6"]], 0)
})

test_that("fitGeneModel can include genome abundance offsets and keep models", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitGeneModel(
        x,
        variable = "condition",
        genomeOffset = TRUE,
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(
        fitInfo(fit)$specification,
        "library_plus_genome_abundance"
    )
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeAssay, "dna_genome_counts")
    expect_identical(fitInfo(fit)$offsetPseudocount, 1)
    expect_identical(sort(names(modelObjects(fit))), sort(rownames(fit)))
    expect_true(all(vapply(
        modelObjects(fit),
        function(x) methods::is(x, "glmmTMB"),
        logical(1)
    )))
})

test_that("fitGeneModel supports numeric sample-level variables", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitGeneModel(x, variable = "pH")

    expect_s4_class(fit, "MTTKFit")
    expect_true(all(as.character(fit$status) == "ok"))
    expect_true(all(as.character(fit$tested_term) == "pH"))
    expect_identical(fitInfo(fit)$variableType, "numeric")
    expect_true(all(is.na(fit$reference_level)))
    expect_true(all(is.na(fit$contrast_level)))
    expect_true(all(as.character(fit$effect_label) == "pH per unit increase"))
})

test_that("fitGeneModel validates modeling inputs", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()

    expect_error(
        fitGeneModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitGeneModel(
            x,
            variable = "condition",
            genomeOffset = TRUE,
            genomeAssay = "missing_assay"
        ),
        "Unknown genome assay name"
    )

    SummarizedExperiment::colData(x)$phase <- factor(c("a", "b", "c", "d"))
    expect_error(
        fitGeneModel(x, variable = "phase"),
        "factor with exactly two levels"
    )
})
