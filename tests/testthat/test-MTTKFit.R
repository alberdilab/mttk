make_test_fit <- function() {
    MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111", "K00368"),
            estimate = c(0.5, -0.3, -1.0),
            p_value = c(0.20, 0.01, 1e-04),
            q_value = c(0.20, 0.03, 1e-04),
            status = c("ok", "ok", "skipped"),
            row.names = c("K03043", "K02111", "K00368")
        ),
        info = list(
            model = "ko_mixed_model",
            effectLabel = "condition: treated vs control"
        ),
        models = list(
            K03043 = list(id = "m1"),
            K02111 = list(id = "m2"),
            K00368 = list(id = "m3")
        )
    )
}

test_that("MTTKFit subsetting prunes stored model objects", {
    fit <- make_test_fit()
    subset_fit <- fit[c("K00368", "K03043"), c("ko_id", "estimate"), drop = FALSE]

    expect_s4_class(subset_fit, "MTTKFit")
    expect_identical(rownames(subset_fit), c("K00368", "K03043"))
    expect_identical(names(modelObjects(subset_fit)), c("K00368", "K03043"))
    expect_identical(modelObjects(subset_fit)[[1]]$id, "m3")
    expect_identical(modelObjects(subset_fit)[[2]]$id, "m1")
})

test_that("koGenomeEffects returns stored KO-by-genome effects and subsets with the fit", {
    fit <- MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111"),
            estimate = c(0.5, -0.3),
            row.names = c("K03043", "K02111")
        ),
        info = list(
            model = "ko_random_slope_model",
            featureIdColumn = "ko_id",
            groupEffectColumn = "genome_id"
        ),
        models = list(
            K03043 = list(id = "m1"),
            K02111 = list(id = "m2")
        )
    )

    metadata_list <- S4Vectors::metadata(fit)
    metadata_list$mttk_fit$groupEffects <- S4Vectors::DataFrame(
        ko_id = c("K03043", "K03043", "K02111"),
        genome_id = c("genome_1", "genome_2", "genome_1"),
        conditional_effect_estimate = c(0.6, 0.4, -0.2),
        row.names = c("K03043::genome_1", "K03043::genome_2", "K02111::genome_1")
    )
    S4Vectors::metadata(fit) <- metadata_list

    effects <- koGenomeEffects(fit)
    expect_s4_class(effects, "DataFrame")
    expect_identical(rownames(effects), c("K03043::genome_1", "K03043::genome_2", "K02111::genome_1"))
    expect_identical(as.character(effects$ko_id), c("K03043", "K03043", "K02111"))

    subset_fit <- fit["K03043", , drop = FALSE]
    subset_effects <- koGenomeEffects(subset_fit)
    expect_identical(rownames(subset_effects), c("K03043::genome_1", "K03043::genome_2"))
    expect_identical(names(modelObjects(subset_fit)), "K03043")
})

test_that("moduleGenomeEffects and pathwayGenomeEffects return stored feature-by-genome effects", {
    module_fit <- MTTKFit(
        results = data.frame(
            module_id = "M00001",
            estimate = 0.5,
            row.names = "M00001"
        ),
        info = list(
            model = "module_random_slope_model",
            featureIdColumn = "module_id",
            groupEffectColumn = "genome_id"
        )
    )
    module_metadata <- S4Vectors::metadata(module_fit)
    module_metadata$mttk_fit$groupEffects <- S4Vectors::DataFrame(
        module_id = "M00001",
        genome_id = "genome_1",
        conditional_effect_estimate = 0.6,
        conditional_effect_std_error = 0.2,
        row.names = "M00001::genome_1"
    )
    S4Vectors::metadata(module_fit) <- module_metadata

    pathway_fit <- MTTKFit(
        results = data.frame(
            pathway_id = "map00010",
            estimate = 0.5,
            row.names = "map00010"
        ),
        info = list(
            model = "pathway_random_slope_model",
            featureIdColumn = "pathway_id",
            groupEffectColumn = "genome_id"
        )
    )
    pathway_metadata <- S4Vectors::metadata(pathway_fit)
    pathway_metadata$mttk_fit$groupEffects <- S4Vectors::DataFrame(
        pathway_id = "map00010",
        genome_id = "genome_1",
        conditional_effect_estimate = 0.6,
        conditional_effect_std_error = 0.2,
        row.names = "map00010::genome_1"
    )
    S4Vectors::metadata(pathway_fit) <- pathway_metadata

    expect_identical(rownames(moduleGenomeEffects(module_fit)), "M00001::genome_1")
    expect_identical(rownames(pathwayGenomeEffects(pathway_fit)), "map00010::genome_1")
})

test_that("fitTable filters and sorts fit results", {
    fit <- make_test_fit()
    tbl <- fitTable(fit, status = "ok", sortBy = "q_value")

    expect_s4_class(tbl, "DataFrame")
    expect_identical(rownames(tbl), c("K02111", "K03043"))
    expect_identical(as.character(tbl$ko_id), c("K02111", "K03043"))

    top_tbl <- fitTable(fit, sortBy = "abs_estimate", n = 1L)
    expect_identical(rownames(top_tbl), "K00368")
})

test_that("coefTable and fitTable(term=) use stored coefficient summaries", {
    fit <- MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111"),
            tested_term = c("conditiontreated", "conditiontreated"),
            estimate = c(0.5, -0.3),
            p_value = c(0.2, 0.01),
            q_value = c(0.2, 0.02),
            status = c("ok", "ok"),
            row.names = c("K03043", "K02111")
        ),
        info = list(
            model = "ko_mixed_model",
            featureIdColumn = "ko_id",
            testedTerm = "conditiontreated"
        ),
        coefficients = data.frame(
            ko_id = c("K03043", "K03043", "K02111", "K02111"),
            model_term = c("conditiontreated", "salinity", "conditiontreated", "salinity"),
            term_variable = c("condition", "salinity", "condition", "salinity"),
            term_variable_type = c("two_level_factor", "numeric", "two_level_factor", "numeric"),
            term_reference_level = c("control", NA, "control", NA),
            term_contrast_level = c("treated", NA, "treated", NA),
            term_effect_label = c(
                "condition: treated vs control",
                "salinity per unit increase",
                "condition: treated vs control",
                "salinity per unit increase"
            ),
            estimate = c(0.5, 0.1, -0.3, -0.4),
            std_error = c(0.2, 0.05, 0.15, 0.08),
            statistic = c(2.5, 2, -2, -5),
            p_value = c(0.2, 0.04, 0.01, 0.001),
            row.names = c(
                "K03043::conditiontreated",
                "K03043::salinity",
                "K02111::conditiontreated",
                "K02111::salinity"
            )
        )
    )

    coefficients <- coefTable(fit, term = "salinity")
    expect_s4_class(coefficients, "DataFrame")
    expect_identical(as.character(coefficients$model_term), c("salinity", "salinity"))
    expect_true("q_value" %in% names(coefficients))

    salinity_table <- fitTable(fit, term = "salinity", sortBy = "p_value", view = "full")
    expect_identical(rownames(salinity_table), c("K02111", "K03043"))
    expect_true(all(as.character(salinity_table$tested_term) == "salinity"))
    expect_equal(as.numeric(salinity_table$estimate), c(-0.4, 0.1))
})

test_that("termFits materializes one canonical fit per stored coefficient term", {
    fit <- MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111"),
            tested_term = c("conditiontreated", "conditiontreated"),
            estimate = c(0.5, -0.3),
            p_value = c(0.2, 0.01),
            q_value = c(0.2, 0.02),
            status = c("ok", "ok"),
            row.names = c("K03043", "K02111")
        ),
        info = list(
            model = "ko_mixed_model",
            featureIdColumn = "ko_id",
            testedTerm = "conditiontreated",
            testedTermInput = "condition",
            effectLabel = "condition: treated vs control"
        ),
        models = list(
            K03043 = list(id = "m1"),
            K02111 = list(id = "m2")
        ),
        coefficients = data.frame(
            ko_id = c("K03043", "K03043", "K02111", "K02111"),
            model_term = c("conditiontreated", "salinity", "conditiontreated", "salinity"),
            term_variable = c("condition", "salinity", "condition", "salinity"),
            term_variable_type = c("two_level_factor", "numeric", "two_level_factor", "numeric"),
            term_reference_level = c("control", NA, "control", NA),
            term_contrast_level = c("treated", NA, "treated", NA),
            term_effect_label = c(
                "condition: treated vs control",
                "salinity per unit increase",
                "condition: treated vs control",
                "salinity per unit increase"
            ),
            estimate = c(0.5, 0.1, -0.3, -0.4),
            std_error = c(0.2, 0.05, 0.15, 0.08),
            statistic = c(2.5, 2, -2, -5),
            p_value = c(0.2, 0.04, 0.01, 0.001),
            row.names = c(
                "K03043::conditiontreated",
                "K03043::salinity",
                "K02111::conditiontreated",
                "K02111::salinity"
            )
        )
    )

    split_fits <- termFits(fit)
    expect_identical(names(split_fits), c("conditiontreated", "salinity"))
    expect_true(all(vapply(split_fits, methods::is, logical(1), class2 = "MTTKFit")))
    expect_identical(fitInfo(split_fits[["salinity"]])$testedTerm, "salinity")
    expect_identical(fitInfo(split_fits[["salinity"]])$effectLabel, "salinity per unit increase")
    expect_true(all(as.character(split_fits[["salinity"]]$tested_term) == "salinity"))
    expect_identical(names(modelObjects(split_fits[["salinity"]])), c("K03043", "K02111"))

    selected <- termFits(fit, terms = "salinity")
    expect_identical(names(selected), "salinity")
})

test_that("termFits only retains random-slope group effects for the original focal term", {
    fit <- MTTKFit(
        results = data.frame(
            ko_id = c("K03043", "K02111"),
            tested_term = c("conditiontreated", "conditiontreated"),
            estimate = c(0.5, -0.3),
            status = c("ok", "ok"),
            row.names = c("K03043", "K02111")
        ),
        info = list(
            model = "ko_random_slope_model",
            featureIdColumn = "ko_id",
            groupEffectColumn = "genome_id",
            testedTerm = "conditiontreated",
            effectLabel = "condition: treated vs control"
        ),
        coefficients = data.frame(
            ko_id = c("K03043", "K03043", "K02111", "K02111"),
            model_term = c("conditiontreated", "salinity", "conditiontreated", "salinity"),
            term_variable = c("condition", "salinity", "condition", "salinity"),
            term_variable_type = c("two_level_factor", "numeric", "two_level_factor", "numeric"),
            term_effect_label = c(
                "condition: treated vs control",
                "salinity per unit increase",
                "condition: treated vs control",
                "salinity per unit increase"
            ),
            estimate = c(0.5, 0.1, -0.3, -0.4),
            std_error = c(0.2, 0.05, 0.15, 0.08),
            p_value = c(0.2, 0.04, 0.01, 0.001),
            row.names = c(
                "K03043::conditiontreated",
                "K03043::salinity",
                "K02111::conditiontreated",
                "K02111::salinity"
            )
        )
    )

    fit_metadata <- S4Vectors::metadata(fit)
    fit_metadata$mttk_fit$groupEffects <- S4Vectors::DataFrame(
        ko_id = c("K03043", "K02111"),
        genome_id = c("genome_1", "genome_1"),
        conditional_effect_estimate = c(0.6, -0.2),
        row.names = c("K03043::genome_1", "K02111::genome_1")
    )
    S4Vectors::metadata(fit) <- fit_metadata

    split_fits <- termFits(fit)
    expect_false(is.null(S4Vectors::metadata(split_fits[["conditiontreated"]])$mttk_fit$groupEffects))
    expect_true(is.null(S4Vectors::metadata(split_fits[["salinity"]])$mttk_fit$groupEffects))
})

test_that("significantResults returns an MTTKFit subset", {
    fit <- make_test_fit()

    sig <- significantResults(fit)
    expect_s4_class(sig, "MTTKFit")
    expect_identical(rownames(sig), "K02111")
    expect_identical(names(modelObjects(sig)), "K02111")

    down <- significantResults(fit, alpha = 1, status = NULL, direction = "down")
    expect_identical(rownames(down), c("K00368", "K02111"))
})

test_that("annotateKOFit joins KO-level annotations from an experiment", {
    fit <- make_test_fit()[c("K03043", "K02111"), , drop = FALSE]
    x <- makeExampleMTTKExperiment()
    annotated <- annotateKOFit(fit, x)

    expect_s4_class(annotated, "MTTKFit")
    expect_true(all(c("module_id", "pathway_id") %in% names(annotated)))
    expect_identical(
        as.character(annotated$module_id),
        c("M00001;M00002", "M00009;M00115")
    )
    expect_identical(
        as.character(annotated$pathway_id),
        c("map00010;map00020", "map00020;map00190")
    )
})

test_that("annotateKOFit can append module descriptions from a KEGG module table", {
    fit <- make_test_fit()["K03043", , drop = FALSE]
    x <- makeExampleMTTKExperiment()
    module_table <- data.frame(
        module_id = c("M00001", "M00002"),
        module_name = c(
            "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate",
            "Glycolysis, core module involving three-carbon compounds"
        ),
        module_class = c(
            "Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism",
            "Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism"
        )
    )

    annotated <- annotateKOFit(fit, x, moduleTable = module_table)

    expect_true(all(c("module_name", "module_class") %in% names(annotated)))
    expect_identical(
        as.character(annotated$module_name),
        paste(
            c(
                "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate",
                "Glycolysis, core module involving three-carbon compounds"
            ),
            collapse = ";"
        )
    )
})

test_that("readKEGGModuleTable standardizes a KEGG module table", {
    tf <- tempfile(fileext = ".tsv")
    writeLines(
        c(
            "module\tdefinition\tname\tclass",
            paste(
                c(
                    "M00001",
                    "K00001 K00002",
                    "Example module",
                    "Pathway modules; Example"
                ),
                collapse = "\t"
            )
        ),
        tf
    )

    module_table <- readKEGGModuleTable(tf)

    expect_s4_class(module_table, "DataFrame")
    expect_identical(
        names(module_table),
        c("module_id", "module_name", "module_class", "module_definition")
    )
    expect_identical(as.character(module_table$module_id), "M00001")
    expect_identical(as.character(module_table$module_name), "Example module")
})

test_that("fetchKEGGModuleTable validates inputs before querying KEGG", {
    expect_error(
        fetchKEGGModuleTable(batchSize = 11),
        "'batchSize' must be between 1 and 10."
    )
    expect_error(
        fetchKEGGModuleTable(sleep = -1),
        "'sleep' must be a single non-negative number."
    )
    expect_error(
        fetchKEGGModuleTable(moduleIds = character()),
        "'moduleIds' does not contain any non-empty module identifiers."
    )
})

test_that("fetchKEGGModuleTable errors cleanly when KEGGREST is unavailable", {
    if (requireNamespace("KEGGREST", quietly = TRUE)) {
        testthat::skip("KEGGREST is installed; skipping the missing-package check.")
    }

    expect_error(
        fetchKEGGModuleTable("M00001"),
        "Package 'KEGGREST' is required"
    )
})

test_that("show method includes the modeled effect when available", {
    fit <- make_test_fit()
    output <- capture.output(show(fit))

    expect_true(any(grepl("^class: MTTKFit$", output)))
    expect_true(any(grepl("^effect: condition: treated vs control$", output)))
})
