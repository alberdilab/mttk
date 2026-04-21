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

test_that("fitTable filters and sorts fit results", {
    fit <- make_test_fit()
    tbl <- fitTable(fit, status = "ok", sortBy = "q_value")

    expect_s4_class(tbl, "DataFrame")
    expect_identical(rownames(tbl), c("K02111", "K03043"))
    expect_identical(as.character(tbl$ko_id), c("K02111", "K03043"))

    top_tbl <- fitTable(fit, sortBy = "abs_estimate", n = 1L)
    expect_identical(rownames(top_tbl), "K00368")
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
