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
    expect_identical(as.character(annotated$module_id), c("M001", "M002"))
    expect_identical(as.character(annotated$pathway_id), c("map00970", "map00910"))
})

test_that("show method includes the modeled effect when available", {
    fit <- make_test_fit()
    output <- capture.output(show(fit))

    expect_true(any(grepl("^class: MTTKFit$", output)))
    expect_true(any(grepl("^effect: condition: treated vs control$", output)))
})
