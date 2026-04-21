test_that("fitGenomeGroupMetaAnalysis summarizes genome effects by genomeData columns", {
    fit <- MTTKFit(
        results = data.frame(
            genome_id = c("genome_1", "genome_2", "genome_3"),
            estimate = c(0.5, 0.3, -0.4),
            std_error = c(0.2, 0.25, 0.3),
            status = c("ok", "ok", "ok"),
            row.names = c("genome_1", "genome_2", "genome_3")
        ),
        info = list(
            model = "genome_model",
            backend = "glmmTMB",
            featureIdColumn = "genome_id",
            variable = "condition",
            variableType = "two_level_factor",
            referenceLevel = "control",
            contrastLevel = "treated",
            effectLabel = "condition: treated vs control",
            formula = "rna_count ~ condition + offset(log(lib_size))"
        )
    )

    x <- makeExampleMTTKExperiment()
    meta_fit <- fitGenomeGroupMetaAnalysis(
        fit,
        x,
        group = "domain",
        method = "fixed",
        minGenomes = 2L
    )

    expect_s4_class(meta_fit, "MTTKFit")
    expect_true(all(c("Bacteria", "Archaea") %in% rownames(meta_fit)))
    expect_identical(fitInfo(meta_fit)$model, "genome_group_meta_analysis")
    expect_identical(fitInfo(meta_fit)$featureIdColumn, "domain")
    expect_identical(fitInfo(meta_fit)$sourceFeatureType, "genome")
    expect_identical(fitInfo(meta_fit)$groupColumn, "domain")
    expect_identical(fitInfo(meta_fit)$metaMethod, "fixed")

    bacteria <- meta_fit["Bacteria", , drop = FALSE]
    expect_identical(as.character(bacteria$domain), "Bacteria")
    expect_identical(as.integer(bacteria$n_genomes_mapped), 2L)
    expect_identical(as.integer(bacteria$n_genomes_tested), 2L)
    expect_identical(as.character(bacteria$member_genome_ids), "genome_1;genome_2")
    expect_identical(as.character(bacteria$status), "ok")
    expect_true(is.finite(as.numeric(bacteria$estimate)))
    expect_identical(as.integer(bacteria$Q_df), 1L)

    archaea <- meta_fit["Archaea", , drop = FALSE]
    expect_identical(as.character(archaea$status), "skipped")
    expect_identical(as.integer(archaea$n_genomes_mapped), 1L)
})

test_that("fitGenomeGroupMetaAnalysis can filter by fit status", {
    fit <- MTTKFit(
        results = data.frame(
            genome_id = c("genome_1", "genome_2", "genome_3"),
            estimate = c(0.5, 0.3, -0.4),
            std_error = c(0.2, 0.25, 0.3),
            status = c("ok", "skipped", "ok"),
            row.names = c("genome_1", "genome_2", "genome_3")
        ),
        info = list(
            model = "genome_model",
            featureIdColumn = "genome_id",
            variable = "condition",
            variableType = "two_level_factor",
            referenceLevel = "control",
            contrastLevel = "treated",
            effectLabel = "condition: treated vs control"
        )
    )

    x <- makeExampleMTTKExperiment()
    meta_fit <- fitGenomeGroupMetaAnalysis(
        fit,
        x,
        group = "domain",
        method = "fixed",
        minGenomes = 1L,
        status = "ok"
    )

    bacteria <- meta_fit["Bacteria", , drop = FALSE]
    expect_identical(as.integer(bacteria$n_genomes_available), 2L)
    expect_identical(as.integer(bacteria$n_genomes_tested), 1L)
    expect_identical(as.character(bacteria$tested_genome_ids), "genome_1")
})

test_that("fitGenomeGroupMetaAnalysis validates fit and grouping inputs", {
    bad_fit <- MTTKFit(
        results = data.frame(
            gene_id = c("gene_1", "gene_2"),
            estimate = c(0.2, -0.1),
            std_error = c(0.1, 0.2),
            row.names = c("gene_1", "gene_2")
        ),
        info = list(
            model = "gene_model",
            featureIdColumn = "gene_id"
        )
    )

    x <- makeExampleMTTKExperiment()

    expect_error(
        fitGenomeGroupMetaAnalysis(bad_fit, x, group = "domain"),
        "genome-level results"
    )
    expect_error(
        fitGenomeGroupMetaAnalysis(
            MTTKFit(
                results = data.frame(
                    genome_id = c("genome_1", "genome_2"),
                    estimate = c(0.2, 0.1),
                    std_error = c(0.1, 0.2),
                    row.names = c("genome_1", "genome_2")
                ),
                info = list(
                    model = "genome_model",
                    featureIdColumn = "genome_id"
                )
            ),
            x,
            group = "missing_group"
        ),
        "Unknown genome grouping column"
    )
})

test_that("fitKOGenomeGroupMetaAnalysis summarizes KO genome effects by genomeData columns", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    meta_fit <- fitKOGenomeGroupMetaAnalysis(
        ko_fit,
        x,
        group = "domain",
        method = "fixed",
        minGenomes = 1L
    )

    expect_s4_class(meta_fit, "MTTKFit")
    expect_identical(fitInfo(meta_fit)$model, "ko_genome_group_meta_analysis")
    expect_identical(fitInfo(meta_fit)$featureIdColumn, "ko_group_id")
    expect_identical(fitInfo(meta_fit)$sourceFitModel, "ko_random_slope_model")
    expect_identical(fitInfo(meta_fit)$groupColumn, "domain")
    expect_identical(fitInfo(meta_fit)$metaMethod, "fixed")
    expect_identical(fitInfo(meta_fit)$effectSource, "conditional_effect_estimate")
    expect_identical(fitInfo(meta_fit)$stdErrorSource, "conditional_effect_std_error")
    expect_true(all(c("ko_id", "domain", "estimate", "n_genomes_tested") %in% names(meta_fit)))
    expect_true(all(grepl("::", rownames(meta_fit), fixed = TRUE)))
    expect_true(any(as.character(meta_fit$domain) == "Bacteria"))
    expect_true(any(as.character(meta_fit$domain) == "Archaea"))
    expect_true(all(as.character(meta_fit$status) == "ok"))
})

test_that("fitKOGenomeGroupMetaAnalysis validates fit and grouping inputs", {
    x <- makeExampleMTTKExperiment()

    bad_fit <- MTTKFit(
        results = data.frame(
            ko_id = "K03043",
            estimate = 0.2,
            std_error = 0.1,
            row.names = "K03043"
        ),
        info = list(
            model = "ko_mixed_model",
            featureIdColumn = "ko_id"
        )
    )

    expect_error(
        fitKOGenomeGroupMetaAnalysis(bad_fit, x, group = "domain"),
        "fitKORandomSlopeModel"
    )
    expect_error(
        fitKOGenomeGroupMetaAnalysis(
            fitKORandomSlopeModel(x, variable = "condition", genomeOffset = FALSE),
            x,
            group = "missing_group"
        ),
        "Unknown genome grouping column"
    )
})

test_that("compareKOGenomeGroups compares KO group responses between groups", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    ko_domain_fit <- fitKOGenomeGroupMetaAnalysis(
        ko_fit,
        x,
        group = "domain",
        method = "fixed",
        minGenomes = 1L
    )

    comparison_fit <- compareKOGenomeGroups(
        ko_domain_fit,
        reference = "Bacteria",
        contrast = "Archaea"
    )

    expect_s4_class(comparison_fit, "MTTKFit")
    expect_identical(fitInfo(comparison_fit)$model, "ko_genome_group_comparison")
    expect_identical(fitInfo(comparison_fit)$groupColumn, "domain")
    expect_identical(fitInfo(comparison_fit)$comparisonMode, "single_contrast")
    expect_identical(fitInfo(comparison_fit)$requestedReference, "Bacteria")
    expect_identical(fitInfo(comparison_fit)$requestedContrast, "Archaea")
    expect_true(all(c(
        "ko_id",
        "reference_domain",
        "contrast_domain",
        "group_comparison",
        "reference_group_estimate",
        "contrast_group_estimate",
        "estimate",
        "std_error"
    ) %in% names(comparison_fit)))
    expect_identical(
        unique(as.character(comparison_fit$reference_domain)),
        "Bacteria"
    )
    expect_identical(
        unique(as.character(comparison_fit$contrast_domain)),
        "Archaea"
    )
    expect_true(any(as.character(comparison_fit$status) == "ok"))
    expect_true(any(as.character(comparison_fit$status) == "skipped"))
    expect_true(all(grepl("::Archaea::Bacteria$", rownames(comparison_fit))))
})

test_that("compareKOGenomeGroups can compare one reference group against all others", {
    fit <- MTTKFit(
        results = data.frame(
            ko_group_id = c(
                "K03043::A", "K03043::B", "K03043::C",
                "K02111::A", "K02111::B", "K02111::C"
            ),
            ko_id = c("K03043", "K03043", "K03043", "K02111", "K02111", "K02111"),
            clade = c("A", "B", "C", "A", "B", "C"),
            estimate = c(0.2, 0.5, -0.1, -0.3, 0.1, -0.2),
            std_error = c(0.1, 0.2, 0.15, 0.12, 0.18, 0.14),
            status = c("ok", "ok", "ok", "ok", "ok", "ok"),
            n_genomes_tested = c(2L, 2L, 2L, 2L, 2L, 2L),
            row.names = c(
                "K03043::A", "K03043::B", "K03043::C",
                "K02111::A", "K02111::B", "K02111::C"
            )
        ),
        info = list(
            model = "ko_genome_group_meta_analysis",
            featureIdColumn = "ko_group_id",
            groupColumn = "clade",
            variable = "condition",
            variableType = "two_level_factor",
            referenceLevel = "control",
            contrastLevel = "treated",
            effectLabel = "condition: treated vs control",
            effectSource = "conditional_effect_estimate",
            stdErrorSource = "conditional_effect_std_error"
        )
    )

    comparison_fit <- compareKOGenomeGroups(
        fit,
        reference = "A"
    )

    expect_s4_class(comparison_fit, "MTTKFit")
    expect_identical(fitInfo(comparison_fit)$comparisonMode, "reference_vs_all")
    expect_identical(fitInfo(comparison_fit)$requestedReference, "A")
    expect_identical(nrow(comparison_fit), 4L)
    expect_true(all(as.character(comparison_fit$reference_clade) == "A"))
    expect_identical(
        sort(unique(as.character(comparison_fit$group_comparison))),
        c("B vs A", "C vs A")
    )
})

test_that("compareKOGenomeGroups validates comparison inputs", {
    fit <- MTTKFit(
        results = data.frame(
            ko_group_id = c("K03043::Bacteria", "K03043::Archaea"),
            ko_id = c("K03043", "K03043"),
            domain = c("Bacteria", "Archaea"),
            estimate = c(0.2, -0.1),
            std_error = c(0.1, 0.2),
            status = c("ok", "ok"),
            n_genomes_tested = c(2L, 1L),
            row.names = c("K03043::Bacteria", "K03043::Archaea")
        ),
        info = list(
            model = "ko_genome_group_meta_analysis",
            featureIdColumn = "ko_group_id",
            groupColumn = "domain",
            variable = "condition",
            variableType = "two_level_factor",
            referenceLevel = "control",
            contrastLevel = "treated",
            effectLabel = "condition: treated vs control"
        )
    )

    expect_error(
        compareKOGenomeGroups(fit, reference = "missing"),
        "Unknown reference group"
    )
    expect_error(
        compareKOGenomeGroups(fit, reference = "Bacteria", contrast = "Bacteria"),
        "must differ"
    )
    expect_error(
        compareKOGenomeGroups(
            MTTKFit(
                results = data.frame(
                    ko_id = "K03043",
                    estimate = 0.2,
                    std_error = 0.1,
                    row.names = "K03043"
                ),
                info = list(
                    model = "ko_mixed_model",
                    featureIdColumn = "ko_id"
                )
            ),
            reference = "Bacteria",
            contrast = "Archaea"
        ),
        "fitKOGenomeGroupMetaAnalysis"
    )
})
