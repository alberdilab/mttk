test_that("fitModuleMetaAnalysis summarizes KO-level effects to modules", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKOMixedModel(x, variable = "condition")
    module_fit <- fitModuleMetaAnalysis(ko_fit, x)

    expect_s4_class(module_fit, "MTTKFit")
    expect_identical(
        rownames(module_fit),
        c("M00001", "M00002", "M00009", "M00010", "M00011", "M00115")
    )
    expect_true(all(c(
        "module_id",
        "estimate",
        "std_error",
        "p_value",
        "q_value",
        "n_kos_mapped",
        "n_kos_available",
        "n_kos_tested",
        "member_ko_ids",
        "tested_ko_ids",
        "mapped_weight_sum",
        "tested_weight_sum",
        "Q",
        "I2",
        "tau2",
        "meta_method",
        "membership_mode",
        "status"
    ) %in% names(module_fit)))

    expect_identical(fitInfo(module_fit)$backend, "meta_analysis")
    expect_identical(fitInfo(module_fit)$model, "module_meta_analysis")
    expect_identical(fitInfo(module_fit)$featureType, "module")
    expect_identical(fitInfo(module_fit)$featureIdColumn, "module_id")
    expect_identical(fitInfo(module_fit)$sourceFitModel, "ko_mixed_model")
    expect_identical(fitInfo(module_fit)$membershipLink, "ko_to_module")
    expect_identical(fitInfo(module_fit)$metaMethod, "random")
    expect_identical(fitInfo(module_fit)$membershipMode, "duplicate")
    expect_true(all(as.character(module_fit$status) == "ok"))
    expect_true(all(as.character(module_fit$membership_mode) == "duplicate"))

    expect_identical(as.character(module_fit$member_ko_ids), c(
        "K03043", "K03043", "K02111", "K00368", "K00368", "K02111"
    ))
    expect_true(all(as.integer(module_fit$n_kos_mapped) == 1L))
    expect_true(all(as.integer(module_fit$n_kos_tested) == 1L))

    ko_estimates <- stats::setNames(as.numeric(ko_fit$estimate), rownames(ko_fit))
    module_estimates <- stats::setNames(as.numeric(module_fit$estimate), rownames(module_fit))
    expect_equal(module_estimates[["M00001"]], ko_estimates[["K03043"]])
    expect_equal(module_estimates[["M00002"]], ko_estimates[["K03043"]])
    expect_equal(module_estimates[["M00009"]], ko_estimates[["K02111"]])
    expect_equal(module_estimates[["M00010"]], ko_estimates[["K00368"]])
    expect_equal(module_estimates[["M00115"]], ko_estimates[["K02111"]])
})

test_that("fitModuleMetaAnalysis accepts KO random-slope fits as input", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKORandomSlopeModel(x, variable = "condition", genomeOffset = FALSE)
    module_fit <- fitModuleMetaAnalysis(ko_fit, x)

    expect_s4_class(module_fit, "MTTKFit")
    expect_identical(fitInfo(module_fit)$sourceFitModel, "ko_random_slope_model")
    expect_identical(
        fitInfo(module_fit)$sourceRandomEffects,
        "genome_random_intercept_and_slope"
    )
    expect_identical(fitInfo(module_fit)$sourceGroupEffectColumn, "genome_id")
    expect_true(all(as.character(module_fit$status) == "ok"))
})

test_that("fitPathwayMetaAnalysis can require multiple KO effects per pathway", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKOMixedModel(x, variable = "condition")
    pathway_fit <- fitPathwayMetaAnalysis(
        ko_fit,
        x,
        method = "fixed",
        minKOs = 2L
    )

    expect_s4_class(pathway_fit, "MTTKFit")
    expect_identical(
        rownames(pathway_fit),
        c("map00010", "map00020", "map00190", "map00680", "map00910")
    )
    expect_identical(fitInfo(pathway_fit)$backend, "meta_analysis")
    expect_identical(fitInfo(pathway_fit)$model, "pathway_meta_analysis")
    expect_identical(fitInfo(pathway_fit)$featureType, "pathway")
    expect_identical(fitInfo(pathway_fit)$featureIdColumn, "pathway_id")
    expect_identical(fitInfo(pathway_fit)$membershipLink, "ko_to_pathway")
    expect_identical(fitInfo(pathway_fit)$metaMethod, "fixed")
    expect_identical(fitInfo(pathway_fit)$membershipMode, "duplicate")
    expect_identical(fitInfo(pathway_fit)$minKOs, 2L)

    expect_identical(
        as.character(pathway_fit$status),
        c("skipped", "ok", "skipped", "skipped", "skipped")
    )
    expect_identical(as.integer(pathway_fit$n_kos_mapped), c(1L, 2L, 1L, 1L, 1L))
    expect_identical(as.integer(pathway_fit$n_kos_tested), c(1L, 2L, 1L, 1L, 1L))
    expect_identical(as.character(pathway_fit$member_ko_ids)[2L], "K03043;K02111")
    expect_identical(as.character(pathway_fit$tested_ko_ids)[2L], "K03043;K02111")
    expect_identical(as.integer(pathway_fit$Q_df)[2L], 1L)
    expect_true(is.finite(as.numeric(pathway_fit$Q)[2L]))
    expect_true(is.finite(as.numeric(pathway_fit$estimate)[2L]))
    expect_true(is.finite(as.numeric(pathway_fit$q_value)[2L]))
    expect_true(is.na(as.numeric(pathway_fit$q_value)[1L]))
})

test_that("fitPathwayMetaAnalysis can split overlapping KO memberships", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    ko_fit <- fitKOMixedModel(x, variable = "condition")
    duplicate_fit <- fitPathwayMetaAnalysis(
        ko_fit,
        x,
        method = "fixed",
        minKOs = 2L,
        membershipMode = "duplicate"
    )
    split_fit <- fitPathwayMetaAnalysis(
        ko_fit,
        x,
        method = "fixed",
        minKOs = 2L,
        membershipMode = "split"
    )

    expect_identical(fitInfo(split_fit)$membershipMode, "split")
    expect_identical(as.character(split_fit$membership_mode)[2L], "split")
    expect_equal(as.numeric(split_fit$tested_weight_sum)[2L], 1)
    expect_equal(
        as.numeric(split_fit$estimate)[2L],
        as.numeric(duplicate_fit$estimate)[2L]
    )
    expect_gt(
        as.numeric(split_fit$std_error)[2L],
        as.numeric(duplicate_fit$std_error)[2L]
    )
})

test_that("KO meta-analysis validates the fit and mapping inputs", {
    bad_fit <- MTTKFit(
        results = data.frame(
            gene_id = c("gene_1", "gene_2"),
            estimate = c(0.2, -0.1),
            std_error = c(0.1, 0.2),
            row.names = c("gene_1", "gene_2")
        ),
        info = list(
            featureIdColumn = "gene_id",
            model = "gene_model"
        )
    )
    x <- makeExampleMTTKExperiment()
    ko_like_fit <- MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111"),
            estimate = c(0.2, -0.1),
            std_error = c(0.1, 0.2),
            status = c("ok", "ok"),
            row.names = c("K03043", "K02111")
        ),
        info = list(
            featureIdColumn = "ko_id",
            model = "ko_mixed_model"
        )
    )

    expect_error(
        fitModuleMetaAnalysis(bad_fit, x),
        "must represent KO-level results"
    )

    y <- x
    links(y) <- links(y)[setdiff(names(links(y)), c("ko_to_module", "module_to_pathway"))]
    expect_error(
        fitModuleMetaAnalysis(ko_like_fit, y),
        "Unknown link table: ko_to_module"
    )
    expect_error(
        fitModuleMetaAnalysis(ko_like_fit, x, membershipMode = "exclusive"),
        "No usable KO memberships remained"
    )
})
