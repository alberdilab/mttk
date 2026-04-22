test_that("makeGenomeModelData aggregates gene RNA to genomes by default", {
    x <- makeExampleMTTKExperiment()
    model_data <- makeGenomeModelData(x, variable = "condition")

    expect_s4_class(model_data, "DataFrame")
    expect_identical(nrow(model_data), 12L)
    expect_true(all(c(
        "genome_id",
        "sample_id",
        "rna_count",
        "genome_name",
        "domain",
        "clade",
        "taxonomy",
        "condition",
        "lib_size",
        "genome_abundance",
        "genome_abundance_offset"
    ) %in% names(model_data)))

    model_state <- S4Vectors::metadata(model_data)$mttk_genome_model
    expect_identical(model_state$specification, "library_plus_genome_abundance")
    expect_identical(model_state$libraryOffset, TRUE)
    expect_identical(model_state$genomeOffset, TRUE)
    expect_identical(model_state$responseAssay, "rna_genome_counts")
    expect_identical(model_state$sourceAssay, "rna_gene_counts")
    expect_identical(model_state$sourceLevel, "aggregated_from_gene_assay")
    expect_identical(model_state$libSizeSource, "colSums(rna_genome_counts)")
    expect_identical(model_state$genomeAssay, "dna_genome_counts")

    one_row <- as.data.frame(
        model_data[
            model_data$genome_id == "genome_1" &
                model_data$sample_id == "sample_1",
            ,
            drop = FALSE
        ]
    )
    expect_identical(one_row$genome_name, "Genome A")
    expect_identical(one_row$domain, "Bacteria")
    expect_equal(one_row$rna_count, 160)
    expect_equal(one_row$lib_size, 405)
    expect_equal(one_row$genome_abundance, dnaGenomeCounts(x)["genome_1", "sample_1"])
    expect_equal(one_row$genome_abundance_offset, one_row$genome_abundance + 1)
})

test_that("genome-level workflows can include a repeated-measures sample block", {
    x <- makeShowcaseMTTKExperiment()
    model_data <- makeGenomeModelData(
        x,
        variable = "condition",
        sampleBlock = "station_id",
        genomeOffset = FALSE
    )

    expect_true("sample_block" %in% names(model_data))
    expect_identical(
        S4Vectors::metadata(model_data)$mttk_genome_model$sampleBlock,
        "station_id"
    )

    skip_if_not_installed("glmmTMB")
    fit <- fitGenomeModel(
        x,
        variable = "condition",
        sampleBlock = "station_id",
        genomeOffset = FALSE
    )

    expect_identical(fitInfo(fit)$sampleBlock, "station_id")
    expect_identical(fitInfo(fit)$randomEffects, "sample_block_random_intercept")
})

test_that("genome-level workflows support formulas with explicit tested terms", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    model_data <- makeGenomeModelData(
        x,
        formula = ~ condition + pH,
        term = "condition",
        genomeOffset = FALSE
    )

    model_state <- S4Vectors::metadata(model_data)$mttk_genome_model
    expect_match(model_state$fixedFormula, "condition")
    expect_match(model_state$fixedFormula, "pH")
    expect_identical(model_state$testedTerm, "conditiontreated")
    expect_identical(model_state$testedVariable, "condition")

    fit <- fitGenomeModel(
        x,
        formula = ~ condition + pH,
        term = "condition",
        genomeOffset = FALSE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_true(all(as.character(fit$tested_term) == "conditiontreated"))
    expect_identical(fitInfo(fit)$testedTerm, "conditiontreated")
    expect_identical(fitInfo(fit)$testedVariable, "condition")
    expect_match(fitInfo(fit)$fixedEffectsFormula, "condition")
    expect_match(fitInfo(fit)$fixedEffectsFormula, "pH")
})

test_that("makeGenomeModelData can use a stored genome-level RNA assay", {
    x <- makeExampleMTTKExperiment()
    aggregated <- aggregateToGenome(x, assays = "rna_gene_counts")
    rnaGenomeCounts(x) <- SummarizedExperiment::assay(
        aggregated,
        "rna_genome_counts",
        withDimnames = TRUE
    )

    model_data <- makeGenomeModelData(
        x,
        variable = "condition",
        assay = "rna_genome_counts",
        genomeOffset = FALSE
    )

    model_state <- S4Vectors::metadata(model_data)$mttk_genome_model
    expect_identical(model_state$specification, "library_only")
    expect_identical(model_state$libraryOffset, TRUE)
    expect_identical(model_state$genomeOffset, FALSE)
    expect_identical(model_state$responseAssay, "rna_genome_counts")
    expect_identical(model_state$sourceAssay, "rna_genome_counts")
    expect_identical(model_state$sourceLevel, "genome_assay")
})

test_that("fitGenomeModel fits one model per genome", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitGenomeModel(x, variable = "condition", genomeOffset = FALSE)

    expect_s4_class(fit, "MTTKFit")
    expect_identical(rownames(fit), c("genome_1", "genome_2", "genome_3"))
    expect_true(all(c(
        "genome_id",
        "genome_name",
        "domain",
        "clade",
        "taxonomy",
        "effect_label",
        "estimate",
        "q_value",
        "status"
    ) %in% names(fit)))
    expect_identical(fitInfo(fit)$backend, "glmmTMB")
    expect_identical(fitInfo(fit)$model, "genome_model")
    expect_identical(fitInfo(fit)$featureType, "genome")
    expect_identical(fitInfo(fit)$featureIdColumn, "genome_id")
    expect_identical(fitInfo(fit)$specification, "library_only")
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, FALSE)
    expect_identical(fitInfo(fit)$responseAssay, "rna_genome_counts")
    expect_identical(fitInfo(fit)$sourceAssay, "rna_gene_counts")
    expect_identical(fitInfo(fit)$sourceLevel, "aggregated_from_gene_assay")
    expect_true(all(as.character(fit$status) == "ok"))

    estimates <- stats::setNames(as.numeric(fit$estimate), rownames(fit))
    expect_gt(estimates[["genome_1"]], 0)
    expect_lt(estimates[["genome_3"]], 0)
})

test_that("fitGenomeModel can use stored genome RNA and genome abundance offsets", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    aggregated <- aggregateToGenome(x, assays = "rna_gene_counts")
    rnaGenomeCounts(x) <- SummarizedExperiment::assay(
        aggregated,
        "rna_genome_counts",
        withDimnames = TRUE
    )

    fit <- fitGenomeModel(
        x,
        variable = "condition",
        assay = "rna_genome_counts",
        genomeOffset = TRUE,
        keepFits = TRUE
    )

    expect_identical(fitInfo(fit)$specification, "library_plus_genome_abundance")
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, TRUE)
    expect_identical(fitInfo(fit)$sourceLevel, "genome_assay")
    expect_identical(fitInfo(fit)$genomeAssay, "dna_genome_counts")
    expect_identical(sort(names(modelObjects(fit))), sort(rownames(fit)))
})

test_that("fitGenomeModel validates modeling inputs", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()

    expect_error(
        fitGenomeModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitGenomeModel(
            x,
            variable = "condition",
            assay = "missing_assay"
        ),
        "Unknown assay name"
    )
    expect_error(
        fitGenomeModel(
            x,
            variable = "condition",
            genomeOffset = TRUE,
            genomeAssay = "missing_assay"
        ),
        "Unknown genome assay name"
    )

    SummarizedExperiment::colData(x)$phase <- factor(c("a", "b", "c", "d"))
    expect_error(
        fitGenomeModel(x, variable = "phase"),
        "produces 3 contrasts"
    )
})
