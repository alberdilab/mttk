test_that("fitKOMixedModel fits the library-only KO mixed model", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKOMixedModel(x, variable = "condition")

    expect_s4_class(fit, "MTTKFit")
    expect_identical(rownames(fit), c("K03043", "K02111", "K00368"))
    expect_true(all(c(
        "ko_id",
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
