test_that("fitGenomePhylogeneticSignal summarizes overall tree signal of genome responses", {
    x <- makeShowcaseMTTKExperiment()
    genome_fit <- fitGenomeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    signal_fit <- fitGenomePhylogeneticSignal(
        genome_fit,
        x,
        nPerm = 19
    )

    expect_s4_class(signal_fit, "MTTKFit")
    expect_identical(rownames(signal_fit), "genome_response_signal")
    expect_identical(signal_fit$analysis_id, "genome_response_signal")
    expect_identical(fitInfo(signal_fit)$model, "genome_phylogenetic_signal")
    expect_true(signal_fit$status %in% c("ok", "skipped"))
    expect_true(is.numeric(signal_fit$signal_statistic))
})

test_that("scanGenomeClades identifies rooted subtree contrasts on genome fits", {
    x <- makeShowcaseMTTKExperiment()
    genome_fit <- fitGenomeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    clade_fit <- scanGenomeClades(
        genome_fit,
        x,
        nPerm = 19,
        minTips = 2
    )

    expect_s4_class(clade_fit, "MTTKFit")
    expect_identical(fitInfo(clade_fit)$model, "genome_clade_scan")
    expect_true(all(c("clade_id", "member_genome_ids", "estimate") %in% names(clade_fit)))
    expect_true(any(clade_fit$member_genome_ids == "genome_1;genome_2"))
    expect_true(any(clade_fit$member_genome_ids == "genome_3;genome_4"))
})

test_that("fitKOPhylogeneticSignal summarizes overall tree signal of KO genome effects", {
    x <- makeShowcaseMTTKExperiment()
    ko_fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    signal_fit <- fitKOPhylogeneticSignal(
        ko_fit,
        x,
        nPerm = 19,
        minGenomes = 3
    )

    expect_s4_class(signal_fit, "MTTKFit")
    expect_identical(fitInfo(signal_fit)$model, "ko_phylogenetic_signal")
    expect_true(all(rownames(signal_fit) == signal_fit$ko_id))
    expect_true(all(c("signal_statistic", "p_value", "q_value") %in% names(signal_fit)))
})

test_that("scanKOClades identifies KO-specific responding clades", {
    x <- makeShowcaseMTTKExperiment()
    ko_fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    clade_fit <- scanKOClades(
        ko_fit,
        x,
        nPerm = 19,
        minTips = 2,
        minGenomes = 3
    )

    expect_s4_class(clade_fit, "MTTKFit")
    expect_identical(fitInfo(clade_fit)$model, "ko_clade_scan")
    expect_true(all(c("ko_clade_id", "ko_id", "clade_id", "member_genome_ids") %in% names(clade_fit)))
    expect_true(all(rownames(clade_fit) == clade_fit$ko_clade_id))
    expect_true(any(clade_fit$member_genome_ids == "genome_1;genome_2"))
})

test_that("module and pathway phylogenetic signal workflows summarize random-slope effects", {
    x <- makeShowcaseMTTKExperiment()
    module_fit <- fitModuleRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    pathway_fit <- fitPathwayRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    module_signal <- fitModulePhylogeneticSignal(
        module_fit,
        x,
        nPerm = 19,
        minGenomes = 3
    )
    pathway_signal <- fitPathwayPhylogeneticSignal(
        pathway_fit,
        x,
        nPerm = 19,
        minGenomes = 3
    )

    expect_s4_class(module_signal, "MTTKFit")
    expect_s4_class(pathway_signal, "MTTKFit")
    expect_identical(fitInfo(module_signal)$model, "module_phylogenetic_signal")
    expect_identical(fitInfo(pathway_signal)$model, "pathway_phylogenetic_signal")
    expect_true(all(rownames(module_signal) == module_signal$module_id))
    expect_true(all(rownames(pathway_signal) == pathway_signal$pathway_id))
    expect_true(all(c("signal_statistic", "p_value", "q_value") %in% names(module_signal)))
    expect_true(all(c("signal_statistic", "p_value", "q_value") %in% names(pathway_signal)))
})

test_that("module and pathway clade scans identify feature-specific responding subtrees", {
    x <- makeShowcaseMTTKExperiment()
    module_fit <- fitModuleRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    pathway_fit <- fitPathwayRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    module_clades <- scanModuleClades(
        module_fit,
        x,
        nPerm = 19,
        minTips = 2,
        minGenomes = 3
    )
    pathway_clades <- scanPathwayClades(
        pathway_fit,
        x,
        nPerm = 19,
        minTips = 2,
        minGenomes = 3
    )

    expect_s4_class(module_clades, "MTTKFit")
    expect_s4_class(pathway_clades, "MTTKFit")
    expect_identical(fitInfo(module_clades)$model, "module_clade_scan")
    expect_identical(fitInfo(pathway_clades)$model, "pathway_clade_scan")
    expect_true(all(c("module_clade_id", "module_id", "clade_id", "member_genome_ids") %in% names(module_clades)))
    expect_true(all(c("pathway_clade_id", "pathway_id", "clade_id", "member_genome_ids") %in% names(pathway_clades)))
    expect_true(any(module_clades$member_genome_ids == "genome_1;genome_2"))
    expect_true(any(pathway_clades$member_genome_ids == "genome_1;genome_2"))
})

test_that("phylogenetic GLS summarizes genome and feature responses under Brownian covariance", {
    x <- makeShowcaseMTTKExperiment()
    genome_fit <- fitGenomeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    ko_fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    module_fit <- fitModuleRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )
    pathway_fit <- fitPathwayRandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE
    )

    genome_gls <- fitGenomePhylogeneticGLS(
        genome_fit,
        x,
        keepFits = TRUE
    )
    ko_gls <- fitKOPhylogeneticGLS(ko_fit, x)
    module_gls <- fitModulePhylogeneticGLS(module_fit, x)
    pathway_gls <- fitPathwayPhylogeneticGLS(pathway_fit, x)

    expect_s4_class(genome_gls, "MTTKFit")
    expect_s4_class(ko_gls, "MTTKFit")
    expect_s4_class(module_gls, "MTTKFit")
    expect_s4_class(pathway_gls, "MTTKFit")
    expect_identical(fitInfo(genome_gls)$model, "genome_phylogenetic_gls")
    expect_identical(fitInfo(ko_gls)$model, "ko_phylogenetic_gls")
    expect_identical(fitInfo(module_gls)$model, "module_phylogenetic_gls")
    expect_identical(fitInfo(pathway_gls)$model, "pathway_phylogenetic_gls")
    expect_identical(rownames(genome_gls), "genome_phylogenetic_gls")
    expect_true(all(c("estimate", "std_error", "statistic", "p_value", "q_value") %in% names(genome_gls)))
    expect_true(all(c("estimate", "std_error", "statistic", "p_value", "q_value") %in% names(ko_gls)))
    expect_true(all(c("estimate", "std_error", "statistic", "p_value", "q_value") %in% names(module_gls)))
    expect_true(all(c("estimate", "std_error", "statistic", "p_value", "q_value") %in% names(pathway_gls)))
    expect_true(length(modelObjects(genome_gls)) == 1L)
    expect_true(inherits(modelObjects(genome_gls)[[1]], "gls"))
})
