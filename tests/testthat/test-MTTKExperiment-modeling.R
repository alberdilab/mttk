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
        genomeOffset = TRUE
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

test_that("aggregateToModuleGenome returns module/genome counts with many-to-many mappings", {
    x <- makeExampleMTTKExperiment()
    aggregated <- aggregateToModuleGenome(x)

    expect_s4_class(aggregated, "SummarizedExperiment")
    expect_identical(nrow(aggregated), 12L)
    expect_true(all(c("module_id", "genome_id", "n_genes_pair", "n_links_pair") %in%
        names(SummarizedExperiment::rowData(aggregated))))
    expect_true(all(c(
        "M00001::genome_1",
        "M00002::genome_1",
        "M00009::genome_3",
        "M00011::genome_2"
    ) %in% rownames(aggregated)))
    expect_equal(
        as.numeric(SummarizedExperiment::assay(aggregated, "rna_gene_counts")["M00001::genome_1", ]),
        as.numeric(rnaGeneCounts(x)["gene_1", ])
    )
    expect_equal(
        as.numeric(SummarizedExperiment::assay(aggregated, "rna_gene_counts")["M00009::genome_3", ]),
        as.numeric(rnaGeneCounts(x)["gene_5", ])
    )
    expect_identical(
        S4Vectors::metadata(aggregated)$mttk_function_aggregation$membershipMode,
        "duplicate"
    )
})

test_that("aggregateToModuleGenome can restrict to exclusive memberships", {
    x <- makeExampleMTTKExperiment()
    link_list <- links(x)
    link_list$ko_to_module <- S4Vectors::DataFrame(
        ko_id = c("K03043", "K02111", "K00368"),
        module_id = c("M00001", "M00009", "M00010")
    )
    link_list$module_to_pathway <- link_list$module_to_pathway[
        as.character(link_list$module_to_pathway$module_id) %in% c("M00001", "M00009", "M00010"),
        ,
        drop = FALSE
    ]
    links(x) <- link_list

    aggregated <- aggregateToModuleGenome(x, membershipMode = "exclusive")

    expect_s4_class(aggregated, "SummarizedExperiment")
    expect_identical(
        rownames(aggregated),
        c(
            "M00001::genome_1",
            "M00009::genome_1",
            "M00001::genome_2",
            "M00010::genome_2",
            "M00009::genome_3",
            "M00010::genome_3"
        )
    )
    expect_identical(
        S4Vectors::metadata(aggregated)$mttk_function_aggregation$membershipMode,
        "exclusive"
    )
})

test_that("aggregateToPathwayGenome can combine multiple genes into one pathway/genome pair", {
    x <- makeExampleMTTKExperiment()
    aggregated <- aggregateToPathwayGenome(x)

    expect_s4_class(aggregated, "SummarizedExperiment")
    expect_identical(nrow(aggregated), 11L)
    expect_true(all(c("pathway_id", "genome_id", "n_genes_pair", "n_links_pair") %in%
        names(SummarizedExperiment::rowData(aggregated))))

    pathway_row <- SummarizedExperiment::rowData(aggregated)["map00020::genome_1", , drop = FALSE]
    expect_identical(as.character(pathway_row$pathway_id), "map00020")
    expect_identical(as.character(pathway_row$genome_id), "genome_1")
    expect_identical(as.integer(pathway_row$n_genes_pair), 2L)
    expect_identical(as.integer(pathway_row$n_links_pair), 2L)
    expect_equal(
        as.numeric(SummarizedExperiment::assay(aggregated, "rna_gene_counts")["map00020::genome_1", ]),
        as.numeric(rnaGeneCounts(x)["gene_1", ] + rnaGeneCounts(x)["gene_2", ])
    )
})

test_that("makeModuleMixedModelData aligns module/genome observations and offsets", {
    x <- makeExampleMTTKExperiment()
    model_data <- makeModuleMixedModelData(
        x,
        variable = "condition",
        genomeOffset = TRUE
    )

    expect_s4_class(model_data, "DataFrame")
    expect_identical(nrow(model_data), 48L)
    expect_true(all(c(
        "module_genome_id",
        "sample_id",
        "rna_count",
        "module_id",
        "genome_id",
        "condition",
        "lib_size",
        "genome_abundance",
        "genome_abundance_offset"
    ) %in% names(model_data)))

    model_state <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    expect_identical(model_state$featureLabel, "module")
    expect_identical(model_state$featureIdColumn, "module_id")
    expect_identical(model_state$specification, "library_plus_genome_abundance")
    expect_identical(model_state$libraryOffset, TRUE)
    expect_identical(model_state$genomeOffset, TRUE)
    expect_identical(model_state$path, c("gene_to_ko", "ko_to_module"))
    expect_identical(model_state$membershipMode, "duplicate")
})

test_that("makePathwayMixedModelData auto-enables genome offsets when available", {
    x <- makeExampleMTTKExperiment()
    model_data <- makePathwayMixedModelData(x, variable = "condition")

    expect_s4_class(model_data, "DataFrame")
    expect_identical(nrow(model_data), 44L)
    expect_true(all(c(
        "pathway_genome_id",
        "sample_id",
        "rna_count",
        "pathway_id",
        "genome_id",
        "condition",
        "lib_size",
        "genome_abundance",
        "genome_abundance_offset"
    ) %in% names(model_data)))

    model_state <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    expect_identical(model_state$featureLabel, "pathway")
    expect_identical(model_state$featureIdColumn, "pathway_id")
    expect_identical(model_state$specification, "library_plus_genome_abundance")
    expect_identical(model_state$libraryOffset, TRUE)
    expect_identical(model_state$genomeOffset, TRUE)
    expect_identical(model_state$path, c("gene_to_ko", "ko_to_pathway"))
    expect_identical(model_state$membershipMode, "duplicate")
})

test_that("fitKOMixedModel fits the library-only KO mixed model", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKOMixedModel(x, variable = "condition", genomeOffset = FALSE)

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
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, FALSE)
    expect_identical(fitInfo(fit)$membershipMode, "duplicate")
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

test_that("fitKORandomSlopeModel stores KO-by-genome conditional effects", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitKORandomSlopeModel(
        x,
        variable = "condition",
        genomeOffset = FALSE,
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(rownames(fit), c("K03043", "K02111", "K00368"))
    expect_identical(fitInfo(fit)$model, "ko_random_slope_model")
    expect_identical(fitInfo(fit)$featureType, "KO")
    expect_identical(fitInfo(fit)$featureIdColumn, "ko_id")
    expect_identical(fitInfo(fit)$specification, "library_only")
    expect_identical(fitInfo(fit)$randomEffects, "genome_random_intercept_and_slope")
    expect_true(all(c(
        "genome_intercept_sd",
        "genome_slope_sd",
        "genome_intercept_slope_cor"
    ) %in% names(fit)))

    effects <- koGenomeEffects(fit)
    expect_s4_class(effects, "DataFrame")
    expect_true(all(c(
        "ko_id",
        "genome_id",
        "conditional_intercept_estimate",
        "conditional_intercept_std_error",
        "conditional_effect_estimate",
        "conditional_effect_std_error",
        "random_intercept_std_error",
        "random_effect_std_error",
        "random_intercept_deviation",
        "random_effect_deviation"
    ) %in% names(effects)))
    expect_true(all(grepl("::", rownames(effects), fixed = TRUE)))
    expect_true(all(unique(as.character(effects$ko_id)) %in% rownames(fit)))
})

test_that("fitModuleMixedModel fits one mixed model per module", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitModuleMixedModel(x, variable = "condition")

    expect_s4_class(fit, "MTTKFit")
    expect_identical(
        rownames(fit),
        c("M00001", "M00002", "M00009", "M00115", "M00010", "M00011")
    )
    expect_true(all(c(
        "module_id",
        "effect_label",
        "estimate",
        "q_value",
        "status",
        "n_genes",
        "n_links",
        "n_genomes"
    ) %in% names(fit)))
    expect_identical(fitInfo(fit)$model, "module_mixed_model")
    expect_identical(fitInfo(fit)$featureType, "module")
    expect_identical(fitInfo(fit)$featureIdColumn, "module_id")
    expect_identical(fitInfo(fit)$specification, "library_plus_genome_abundance")
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, TRUE)
    expect_identical(fitInfo(fit)$aggregationPath, c("gene_to_ko", "ko_to_module"))
    expect_identical(fitInfo(fit)$membershipMode, "duplicate")
    expect_true(all(as.character(fit$status) == "ok"))

    estimates <- stats::setNames(as.numeric(fit$estimate), rownames(fit))
    expect_gt(estimates[["M00001"]], 0)
    expect_lt(estimates[["M00010"]], 0)
})

test_that("fitModuleRandomSlopeModel stores module-by-genome conditional effects", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitModuleRandomSlopeModel(
        x,
        variable = "condition",
        membershipMode = "duplicate",
        genomeOffset = FALSE,
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(fitInfo(fit)$model, "module_random_slope_model")
    expect_identical(fitInfo(fit)$featureType, "module")
    expect_identical(fitInfo(fit)$featureIdColumn, "module_id")
    expect_identical(fitInfo(fit)$membershipMode, "duplicate")
    expect_identical(fitInfo(fit)$randomEffects, "genome_random_intercept_and_slope")

    effects <- moduleGenomeEffects(fit)
    expect_s4_class(effects, "DataFrame")
    expect_true(all(c(
        "module_id",
        "genome_id",
        "conditional_effect_estimate",
        "conditional_effect_std_error"
    ) %in% names(effects)))
    expect_true(all(unique(as.character(effects$module_id)) %in% rownames(fit)))
})

test_that("fitPathwayMixedModel fits one mixed model per pathway", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitPathwayMixedModel(
        x,
        variable = "condition",
        genomeOffset = TRUE,
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(
        rownames(fit),
        c("map00010", "map00020", "map00190", "map00910", "map00680")
    )
    expect_true(all(c(
        "pathway_id",
        "effect_label",
        "estimate",
        "q_value",
        "status"
    ) %in% names(fit)))
    expect_identical(fitInfo(fit)$model, "pathway_mixed_model")
    expect_identical(fitInfo(fit)$featureType, "pathway")
    expect_identical(fitInfo(fit)$featureIdColumn, "pathway_id")
    expect_identical(fitInfo(fit)$specification, "library_plus_genome_abundance")
    expect_identical(fitInfo(fit)$libraryOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeOffset, TRUE)
    expect_identical(fitInfo(fit)$genomeAssay, "dna_genome_counts")
    expect_identical(fitInfo(fit)$aggregationPath, c("gene_to_ko", "ko_to_pathway"))
    expect_identical(fitInfo(fit)$membershipMode, "duplicate")
    expect_identical(sort(names(modelObjects(fit))), sort(rownames(fit)))
    expect_true(all(as.character(fit$status) == "ok"))

    estimates <- stats::setNames(as.numeric(fit$estimate), rownames(fit))
    expect_gt(estimates[["map00010"]], 0)
    expect_lt(estimates[["map00910"]], 0)
})

test_that("fitPathwayRandomSlopeModel stores pathway-by-genome conditional effects", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    fit <- fitPathwayRandomSlopeModel(
        x,
        variable = "condition",
        membershipMode = "duplicate",
        genomeOffset = FALSE,
        keepFits = TRUE
    )

    expect_s4_class(fit, "MTTKFit")
    expect_identical(fitInfo(fit)$model, "pathway_random_slope_model")
    expect_identical(fitInfo(fit)$featureType, "pathway")
    expect_identical(fitInfo(fit)$featureIdColumn, "pathway_id")
    expect_identical(fitInfo(fit)$membershipMode, "duplicate")
    expect_identical(fitInfo(fit)$randomEffects, "genome_random_intercept_and_slope")

    effects <- pathwayGenomeEffects(fit)
    expect_s4_class(effects, "DataFrame")
    expect_true(all(c(
        "pathway_id",
        "genome_id",
        "conditional_effect_estimate",
        "conditional_effect_std_error"
    ) %in% names(effects)))
    expect_true(all(unique(as.character(effects$pathway_id)) %in% rownames(fit)))
})

test_that("fitModuleMixedModel can use exclusive memberships when mappings are unique", {
    skip_if_not_installed("glmmTMB")

    x <- makeExampleMTTKExperiment()
    link_list <- links(x)
    link_list$ko_to_module <- S4Vectors::DataFrame(
        ko_id = c("K03043", "K02111", "K00368"),
        module_id = c("M00001", "M00009", "M00010")
    )
    link_list$module_to_pathway <- link_list$module_to_pathway[
        as.character(link_list$module_to_pathway$module_id) %in% c("M00001", "M00009", "M00010"),
        ,
        drop = FALSE
    ]
    links(x) <- link_list

    fit <- fitModuleMixedModel(x, variable = "condition", membershipMode = "exclusive")

    expect_s4_class(fit, "MTTKFit")
    expect_identical(fitInfo(fit)$membershipMode, "exclusive")
    expect_true(all(as.character(fit$status) == "ok"))
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
            genomeOffset = TRUE,
            genomeAssay = "missing_assay"
        ),
        "Unknown genome assay name"
    )

    SummarizedExperiment::colData(x)$phase <- factor(c("a", "b", "c", "d"))
    expect_error(
        fitKOMixedModel(x, variable = "phase"),
        "factor with exactly two levels"
    )
    expect_error(
        fitModuleMixedModel(
            x,
            variable = "condition",
            genomeOffset = TRUE,
            genomeAssay = "missing_assay"
        ),
        "Unknown genome assay name"
    )
    expect_error(
        fitPathwayMixedModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitKORandomSlopeModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitModuleRandomSlopeModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
    expect_error(
        fitPathwayRandomSlopeModel(x, variable = "missing_variable"),
        "Unknown sample-level variable"
    )
})
