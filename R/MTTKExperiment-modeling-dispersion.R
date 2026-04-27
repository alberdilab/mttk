.fit_glmmtmb_nbinom2_disp <- function(formula, dispformula, data, formula_environment = NULL) {
    formula <- .set_formula_environment(formula = formula, values = formula_environment)

    if (is.null(formula_environment) || length(formula_environment) == 0L) {
        return(
            glmmTMB::glmmTMB(
                formula     = formula,
                dispformula = dispformula,
                data        = data,
                family      = glmmTMB::nbinom2(link = "log")
            )
        )
    }

    fit_environment <- list2env(
        c(list(formula = formula, dispformula = dispformula, data = data),
          formula_environment),
        parent = environment(formula)
    )
    eval(
        quote(
            glmmTMB::glmmTMB(
                formula     = formula,
                dispformula = dispformula,
                data        = data,
                family      = glmmTMB::nbinom2(link = "log")
            )
        ),
        envir = fit_environment
    )
}

.empty_ko_dispersion_row <- function(ko_id, ko_summary, effect_label) {
    S4Vectors::DataFrame(
        ko_id                  = ko_id,
        effect_label           = effect_label,
        estimate               = NA_real_,
        std_error              = NA_real_,
        statistic              = NA_real_,
        p_value                = NA_real_,
        q_value                = NA_real_,
        n_genes                = as.integer(ko_summary$n_genes[[1L]]),
        n_genomes              = as.integer(ko_summary$n_genomes[[1L]]),
        n_observations         = NA_integer_,
        n_nonzero_observations = NA_integer_,
        AIC                    = NA_real_,
        BIC                    = NA_real_,
        logLik                 = NA_real_,
        pd_hess                = NA,
        optimizer_convergence  = NA_integer_,
        warning_message        = NA_character_,
        error_message          = NA_character_,
        status                 = NA_character_,
        row.names              = ko_id
    )
}

.fit_one_ko_dispersion_model <- function(data, ko_id, ko_summary,
                                          model_spec, mean_formula,
                                          dispformula, formula_environment,
                                          keep_fits) {
    effect_label <- paste0("dispersion: ", model_spec$effectLabel)
    row <- .empty_ko_dispersion_row(ko_id = ko_id, ko_summary = ko_summary,
                                    effect_label = effect_label)
    row$n_observations         <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0L)

    if (nrow(data) == 0L) {
        row$status        <- "skipped"
        row$error_message <- "No observations were available for this KO."
        return(list(result = row, model = NULL))
    }
    if (length(unique(as.character(data$genome_id))) < 2L) {
        row$status        <- "skipped"
        row$error_message <- "At least two genomes are required to estimate the random effect."
        return(list(result = row, model = NULL))
    }
    if (all(data$rna_count == 0L)) {
        row$status        <- "skipped"
        row$error_message <- "All KO-level counts were zero."
        return(list(result = row, model = NULL))
    }

    warning_messages <- character()
    fit <- tryCatch(
        withCallingHandlers(
            .fit_glmmtmb_nbinom2_disp(
                formula            = mean_formula,
                dispformula        = dispformula,
                data               = data,
                formula_environment = formula_environment
            ),
            warning = function(w) {
                warning_messages <<- c(warning_messages, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = identity
    )

    if (inherits(fit, "error")) {
        row$status        <- "error"
        row$error_message <- conditionMessage(fit)
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else NA_character_
        return(list(result = row, model = NULL))
    }

    disp_table <- tryCatch(
        summary(fit)$coefficients$disp,
        error = function(e) NULL
    )
    if (is.null(disp_table) || nrow(disp_table) < 2L) {
        row$status        <- "error"
        row$error_message <- "Dispersion coefficient table could not be extracted or had only an intercept."
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else NA_character_
        return(list(result = row, model = if (keep_fits) fit else NULL))
    }

    tested_row_name <- .match_model_term_name(
        model_spec$testedTerm,
        rownames(disp_table)
    )
    if (is.na(tested_row_name)) {
        non_intercept <- setdiff(rownames(disp_table), "(Intercept)")
        tested_row_name <- if (length(non_intercept) > 0L) non_intercept[[1L]] else NA_character_
    }
    if (is.na(tested_row_name)) {
        row$status        <- "error"
        row$error_message <- "Could not identify the tested dispersion term in the coefficient table."
        return(list(result = row, model = if (keep_fits) fit else NULL))
    }

    statistic_col <- intersect(colnames(disp_table), c("z value", "t value"))[[1L]]
    p_value_col   <- grep("^Pr\\(", colnames(disp_table), value = TRUE)[[1L]]
    tested_row    <- disp_table[tested_row_name, , drop = FALSE]

    row$estimate              <- as.numeric(tested_row[, "Estimate"])
    row$std_error             <- as.numeric(tested_row[, "Std. Error"])
    row$statistic             <- as.numeric(tested_row[, statistic_col])
    row$p_value               <- as.numeric(tested_row[, p_value_col])
    row$AIC                   <- stats::AIC(fit)
    row$BIC                   <- stats::BIC(fit)
    row$logLik                <- as.numeric(stats::logLik(fit))
    row$pd_hess               <- if (!is.null(fit$sdr$pdHess)) isTRUE(fit$sdr$pdHess) else NA
    row$optimizer_convergence <- if (!is.null(fit$fit$convergence)) as.integer(fit$fit$convergence) else NA_integer_
    row$warning_message       <- if (length(warning_messages) > 0L) {
        paste(unique(warning_messages), collapse = " | ")
    } else NA_character_
    row$error_message <- NA_character_
    row$status        <- "ok"

    list(result = row, model = if (keep_fits) fit else NULL)
}

#' Test Whether KO Expression Variability Differs Across Conditions
#'
#' `fitKODispersion()` fits one negative-binomial mixed model per KO with a
#' variable dispersion term (`dispformula = ~ variable`), testing whether the
#' overdispersion parameter — and therefore the sample-to-sample variability of
#' expression — differs across conditions or along a continuous gradient.
#'
#' The mean structure mirrors [fitKOMixedModel()]: RNA counts are aggregated to
#' KO-within-genome observations, a genome random intercept is included, and
#' optional library-size and genome-abundance offsets are applied. The
#' dispersion sub-model adds a separate linear predictor for `log(phi)` where
#' `phi` is the NB2 size parameter.
#'
#' **Interpreting the estimate**: the reported estimate is the log-ratio of the
#' dispersion parameter between the contrast and reference level. A positive
#' value means higher `phi` in the contrast group — i.e. more Poisson-like,
#' *less* overdispersed, *more* consistent expression. A negative value means
#' *more* overdispersed, *more* variable expression in the contrast group.
#'
#' @inheritParams fitKOMixedModel
#'
#' @return An `MTTKFit` with one row per KO. The `estimate` column is the
#'   log-ratio of the dispersion parameter (`log(phi_contrast / phi_reference)`).
#'
#' @seealso [fitKOMixedModel()] for mean-level KO analysis,
#'   [varianceDecomposition()] for decomposing variance into genome and
#'   residual components from a fitted model.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     disp_fit <- fitKODispersion(x, variable = "condition")
#'     fitTable(disp_fit)
#' }
#'
#' @export
fitKODispersion <- function(
    x,
    variable        = NULL,
    formula         = NULL,
    term            = NULL,
    assay           = "rna_gene_counts",
    libraryOffset   = TRUE,
    genomeOffset    = NULL,
    libSize         = NULL,
    sampleBlock     = NULL,
    genomeAssay     = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL,
    keepFits        = FALSE,
    BPPARAM         = BiocParallel::SerialParam()
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use fitKODispersion().",
            call. = FALSE
        )
    }
    if (!is.logical(keepFits) || length(keepFits) != 1L || is.na(keepFits)) {
        stop("'keepFits' must be TRUE or FALSE.", call. = FALSE)
    }

    model_spec <- .normalize_model_spec(
        x               = x,
        variable        = variable,
        formula         = formula,
        term            = term,
        referenceLevels = referenceLevels,
        allow_formula   = TRUE,
        allow_random_slope = FALSE
    )
    .require_model_spec_term(model_spec)

    model_data <- .make_feature_mixed_model_data(
        x              = x,
        model_spec     = model_spec,
        assay          = assay,
        libraryOffset  = libraryOffset,
        genomeOffset   = genomeOffset,
        membershipMode = "duplicate",
        libSize        = libSize,
        sampleBlock    = sampleBlock,
        genomeAssay    = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        path           = "gene_to_ko",
        feature_label  = "KO"
    )

    observations   <- as.data.frame(model_data)
    model_state    <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    specification  <- model_state$specification
    feature_summary <- model_state$featureSummary
    feature_ids    <- rownames(feature_summary)

    random_term  <- "(1 | genome_id)"
    has_block    <- "sample_block" %in% names(observations)
    block_term   <- if (has_block) " + (1 | sample_block)" else ""
    mean_formula <- stats::as.formula(
        paste0("rna_count ~ ", specification, " + offset(log_offset) + ",
               random_term, block_term)
    )
    disp_variable <- model_spec$testedVariable
    if (is.null(disp_variable) || is.na(disp_variable) || disp_variable == "") {
        stop("Cannot determine the dispersion variable from the model specification.",
             call. = FALSE)
    }
    dispformula <- stats::as.formula(paste0("~ `", disp_variable, "`"))

    fitted_rows <- BiocParallel::bplapply(feature_ids, function(ko_id) {
        data_ko <- observations[
            as.character(observations$ko_id) == ko_id, , drop = FALSE
        ]
        .fit_one_ko_dispersion_model(
            data               = data_ko,
            ko_id              = ko_id,
            ko_summary         = feature_summary[ko_id, , drop = FALSE],
            model_spec         = model_spec,
            mean_formula       = mean_formula,
            dispformula        = dispformula,
            formula_environment = NULL,
            keep_fits          = keepFits
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(rbind, lapply(fitted_rows, function(r) r$result))
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(lapply(fitted_rows, function(r) r$model), feature_ids)
    } else {
        list()
    }

    info <- list(
        backend           = "glmmTMB",
        model             = "ko_dispersion_model",
        featureType       = "KO",
        featureIdColumn   = "ko_id",
        variable          = model_state$variable,
        variableType      = model_state$variableType,
        referenceLevel    = model_state$referenceLevel,
        contrastLevel     = model_state$contrastLevel,
        effectLabel       = model_state$effectLabel,
        testedTermInput   = model_state$testedTermInput,
        testedTerm        = model_state$testedTerm,
        testedVariable    = model_state$testedVariable,
        dispersionFormula = paste(deparse(dispformula), collapse = " "),
        family            = "nbinom2",
        libraryOffset     = model_state$libraryOffset,
        genomeOffset      = model_state$genomeOffset,
        sampleBlock       = model_state$sampleBlock,
        n_features        = nrow(results),
        n_ok              = sum(results$status == "ok"),
        n_skipped         = sum(results$status == "skipped"),
        n_error           = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info    = info,
        models  = stored_models,
        coefficients = S4Vectors::DataFrame()
    )
}
