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
