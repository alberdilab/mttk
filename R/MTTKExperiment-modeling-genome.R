.prepare_genome_offset_inputs <- function(
    x,
    response_counts,
    response_assay,
    libraryOffset,
    genomeOffset,
    libSize,
    genomeAssay,
    offsetPseudocount
) {
    if (!is.numeric(offsetPseudocount) ||
        length(offsetPseudocount) != 1L ||
        is.na(offsetPseudocount) ||
        offsetPseudocount < 0) {
        stop("'offsetPseudocount' must be a single non-negative numeric value.", call. = FALSE)
    }

    offset_flags <- .normalize_offset_flags(
        x = x,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        genomeAssay = genomeAssay
    )

    lib_size <- if (offset_flags$libraryOffset) {
        .normalize_genome_model_lib_size(
            x = x,
            response_counts = response_counts,
            response_assay = response_assay,
            libSize = libSize
        )
    } else {
        list(
            values = stats::setNames(rep(NA_real_, ncol(response_counts)), colnames(response_counts)),
            source = NA_character_
        )
    }

    normalized_genome_assay <- if (offset_flags$genomeOffset) {
        .normalize_genome_analysis_assay(x, assay = genomeAssay)
    } else {
        NA_character_
    }

    list(
        libraryOffset = offset_flags$libraryOffset,
        genomeOffset = offset_flags$genomeOffset,
        specification = offset_flags$specification,
        libSize = lib_size,
        genomeAssay = normalized_genome_assay,
        offsetPseudocount = if (offset_flags$genomeOffset) offsetPseudocount else NA_real_
    )
}

.build_genome_model_observations <- function(
    x,
    model_spec,
    specification,
    response_counts,
    genome_summary,
    response_assay,
    response_source_assay,
    response_source_level,
    lib_size,
    block_spec,
    genome_assay,
    offset_pseudocount
) {
    sample_ids <- colnames(response_counts)
    genome_ids <- rownames(response_counts)

    if (!("genome_id" %in% names(genome_summary))) {
        genome_summary$genome_id <- genome_ids
    }
    rownames(genome_summary) <- genome_ids

    obs <- expand.grid(
        genome_id = genome_ids,
        sample_id = sample_ids,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    obs$rna_count <- as.numeric(as.vector(response_counts))

    matched_genomes <- match(obs$genome_id, genome_ids)
    summary_columns <- setdiff(names(genome_summary), "genome_id")
    for (column_name in summary_columns) {
        obs[[column_name]] <- S4Vectors::decode(genome_summary[[column_name]])[matched_genomes]
    }

    predictor_frame <- model_spec$sampleData[obs$sample_id, , drop = FALSE]
    predictor_names <- names(predictor_frame)
    for (column_name in predictor_names) {
        values <- predictor_frame[[column_name]]
        if (is.factor(values)) {
            obs[[column_name]] <- factor(values, levels = levels(values))
        } else {
            obs[[column_name]] <- values
        }
    }

    if (!is.null(block_spec$values)) {
        obs$sample_block <- factor(
            as.character(block_spec$values[obs$sample_id]),
            levels = levels(block_spec$values)
        )
    }

    if (.specification_uses_library_offset(specification)) {
        obs$lib_size <- as.numeric(lib_size$values[obs$sample_id])
    }

    if (.specification_uses_genome_offset(specification)) {
        genome_mat <- SummarizedExperiment::assay(
            genomeExperiment(x),
            genome_assay,
            withDimnames = TRUE
        )
        matched_offset_genomes <- match(obs$genome_id, rownames(genome_mat))
        matched_samples <- match(obs$sample_id, colnames(genome_mat))

        if (anyNA(matched_offset_genomes) || anyNA(matched_samples)) {
            stop(
                "Every genome/sample observation must be present in the selected genome assay.",
                call. = FALSE
            )
        }

        obs$genome_abundance <- as.numeric(genome_mat[cbind(matched_offset_genomes, matched_samples)])
        obs$genome_abundance_offset <- obs$genome_abundance + offset_pseudocount
    }

    out <- S4Vectors::DataFrame(obs, check.names = FALSE)
    S4Vectors::metadata(out)$mttk_genome_model <- list(
        variable = model_spec$variable,
        variableType = model_spec$variableType,
        referenceLevel = model_spec$referenceLevel,
        contrastLevel = model_spec$contrastLevel,
        effectLabel = model_spec$effectLabel,
        fixedFormula = model_spec$fixedFormulaLabel,
        availableTerms = model_spec$availableTerms,
        testedTermInput = model_spec$testedTermInput,
        testedTerm = model_spec$testedTerm,
        testedVariable = model_spec$testedVariable,
        sampleBlock = block_spec$name,
        specification = specification,
        libraryOffset = .specification_uses_library_offset(specification),
        genomeOffset = .specification_uses_genome_offset(specification),
        responseAssay = response_assay,
        sourceAssay = response_source_assay,
        sourceLevel = response_source_level,
        libSizeSource = if (.specification_uses_library_offset(specification)) {
            lib_size$source
        } else {
            NA_character_
        },
        genomeAssay = if (.specification_uses_genome_offset(specification)) {
            genome_assay
        } else {
            NA_character_
        },
        offsetPseudocount = if (.specification_uses_genome_offset(specification)) {
            offset_pseudocount
        } else {
            NA_real_
        },
        genomeSummary = genome_summary
    )

    out
}

.empty_genome_model_row <- function(genome_id, genome_summary, variable_info) {
    columns <- list(genome_id = genome_id)
    extra_columns <- setdiff(names(genome_summary), "genome_id")

    for (column_name in extra_columns) {
        columns[[column_name]] <- S4Vectors::decode(genome_summary[[column_name]])[[1L]]
    }

    do.call(
        S4Vectors::DataFrame,
        c(
            columns,
            list(
                tested_term = NA_character_,
                variable_type = variable_info$type,
                reference_level = variable_info$referenceLevel,
                contrast_level = variable_info$contrastLevel,
                effect_label = variable_info$effectLabel,
                estimate = NA_real_,
                std_error = NA_real_,
                statistic = NA_real_,
                p_value = NA_real_,
                q_value = NA_real_,
                intercept_estimate = NA_real_,
                intercept_std_error = NA_real_,
                n_observations = NA_integer_,
                n_nonzero_observations = NA_integer_,
                AIC = NA_real_,
                BIC = NA_real_,
                logLik = NA_real_,
                pd_hess = NA,
                optimizer_convergence = NA_integer_,
                warning_message = NA_character_,
                error_message = NA_character_,
                status = NA_character_,
                row.names = genome_id
            )
        )
    )
}

.fit_one_genome_model <- function(data, genome_id, genome_summary, variable_info, specification, keep_fits) {
    row <- .empty_genome_model_row(
        genome_id = genome_id,
        genome_summary = genome_summary,
        variable_info = variable_info
    )
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- "No observations were available for the genome."
        return(list(result = row, model = NULL))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- "All genome-level counts were zero."
        return(list(result = row, model = NULL))
    }

    warning_messages <- character()
    formula <- .gene_model_formula(
        model_spec = variable_info$modelSpec,
        specification = specification,
        sample_block = "sample_block" %in% names(data)
    )
    fit <- tryCatch(
        withCallingHandlers(
            glmmTMB::glmmTMB(
                formula = formula,
                data = data,
                family = glmmTMB::nbinom2(link = "log")
            ),
            warning = function(w) {
                warning_messages <<- c(warning_messages, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = identity
    )

    if (inherits(fit, "error")) {
        row$status <- "error"
        row$error_message <- conditionMessage(fit)
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else {
            NA_character_
        }
        return(list(result = row, model = NULL))
    }

    coefficient_table <- summary(fit)$coefficients$cond
    tested_term <- tryCatch(
        .resolve_fitted_tested_term(variable_info$modelSpec, coefficient_table),
        error = identity
    )

    if (inherits(tested_term, "error")) {
        row$status <- "error"
        row$error_message <- conditionMessage(tested_term)
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else {
            NA_character_
        }
        return(list(result = row, model = if (keep_fits) fit else NULL))
    }

    statistic_col <- intersect(colnames(coefficient_table), c("z value", "t value"))[1L]
    p_value_col <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)[1L]
    tested_row <- coefficient_table[tested_term, , drop = FALSE]
    has_intercept <- "(Intercept)" %in% rownames(coefficient_table)
    intercept_row <- if (has_intercept) coefficient_table["(Intercept)", , drop = FALSE] else NULL
    term_info <- .resolved_term_info(variable_info$modelSpec, tested_term = tested_term)

    row$tested_term <- tested_term
    row$variable_type <- term_info$type
    row$reference_level <- term_info$referenceLevel
    row$contrast_level <- term_info$contrastLevel
    row$effect_label <- term_info$effectLabel
    row$estimate <- as.numeric(tested_row[, "Estimate"])
    row$std_error <- as.numeric(tested_row[, "Std. Error"])
    row$statistic <- as.numeric(tested_row[, statistic_col])
    row$p_value <- as.numeric(tested_row[, p_value_col])
    row$intercept_estimate <- if (!is.null(intercept_row)) as.numeric(intercept_row[, "Estimate"]) else NA_real_
    row$intercept_std_error <- if (!is.null(intercept_row)) as.numeric(intercept_row[, "Std. Error"]) else NA_real_
    row$AIC <- stats::AIC(fit)
    row$BIC <- stats::BIC(fit)
    row$logLik <- as.numeric(stats::logLik(fit))
    row$pd_hess <- if (!is.null(fit$sdr$pdHess)) isTRUE(fit$sdr$pdHess) else NA
    row$optimizer_convergence <- if (!is.null(fit$fit$convergence)) {
        as.integer(fit$fit$convergence)
    } else {
        NA_integer_
    }
    row$warning_message <- if (length(warning_messages) > 0L) {
        paste(unique(warning_messages), collapse = " | ")
    } else {
        NA_character_
    }
    row$error_message <- NA_character_
    row$status <- "ok"

    list(
        result = row,
        model = if (keep_fits) fit else NULL
    )
}

#' Build a Genome-Level Model Data Table
#'
#' `makeGenomeModelData()` materializes the long-form observation table used by
#' [fitGenomeModel()]. Each row corresponds to one genome/sample observation,
#' with genome-level RNA counts, sample-level covariates, genome metadata, and
#' an optional genome-abundance offset aligned and ready for model fitting.
#'
#' This workflow is intended for the question of which genomes change across
#' conditions or are associated with a continuous variable.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)` for the
#'   simple one-variable interface. The column must be numeric or a factor with
#'   at least two levels. For factors with three or more levels supply `term` to
#'   select which contrast to test. Supply exactly one of `variable` or `formula`.
#' @param formula Optional one-sided or two-sided fixed-effect formula for the
#'   sample-level covariates, for example `~ condition + pH` or
#'   `rna_count ~ condition + pH`. Offsets are added internally by mttk.
#'   Supply exactly one of `variable` or `formula`.
#' @param term Optional character string naming the model term to test. Required
#'   when `variable` refers to a factor with three or more levels (supply the
#'   contrast name, e.g. `"conditionB"`) or when `formula` defines more than one
#'   tested term.
#' @param assay Genome-level or gene-level assay name used as the RNA response.
#'   Use `NULL` to prefer `"rna_genome_counts"` when present and otherwise
#'   aggregate `"rna_gene_counts"` to genomes.
#' @param libraryOffset Logical; if `TRUE`, include `offset(log(lib_size))`.
#' @param libSize Library-size offset specification. Use `NULL` to compute
#'   `colSums()` of the resolved genome-level response assay, a single
#'   `colData(x)` column name, or a numeric vector with one value per sample.
#' @param sampleBlock Optional sample-level blocking variable from
#'   `colData(x)`. When supplied, mttk adds a random intercept
#'   `(1 | sample_block)` to account for repeated measures or paired samples.
#' @param genomeOffset Logical; if `TRUE`, include
#'   `offset(log(genome_abundance + offsetPseudocount))`. If `NULL` (the
#'   default), mttk uses genome-abundance normalization when `genomeAssay` is
#'   available in `genomeExperiment(x)`.
#' @param genomeAssay Genome-level assay used when `genomeOffset = TRUE`.
#' @param offsetPseudocount Non-negative pseudocount added to genome abundance
#'   before log-offset calculation.
#' @param referenceLevels Optional named list or named character vector setting
#'   the reference levels of factor-like sample covariates before model fitting.
#'
#' @return An `S4Vectors::DataFrame` with one row per genome/sample observation.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' genome_model_data <- makeGenomeModelData(x, variable = "condition")
#' genome_model_data
#'
#' @export
makeGenomeModelData <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = NULL,
    libraryOffset = TRUE,
    genomeOffset = NULL,
    libSize = NULL,
    sampleBlock = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    model_spec <- .normalize_model_spec(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        referenceLevels = referenceLevels,
        allow_formula = TRUE,
        allow_random_slope = FALSE
    )
    block_spec <- .normalize_sample_block_spec(
        x = x,
        sampleBlock = sampleBlock,
        fixed_variables = model_spec$fixedVariables
    )

    response_state <- .resolve_genome_model_response(x, assay = assay)
    offset_state <- .prepare_genome_offset_inputs(
        x = x,
        response_counts = response_state$counts,
        response_assay = response_state$responseAssay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        libSize = libSize,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount
    )

    .build_genome_model_observations(
        x = x,
        model_spec = model_spec,
        specification = offset_state$specification,
        response_counts = response_state$counts,
        genome_summary = response_state$genomeSummary,
        response_assay = response_state$responseAssay,
        response_source_assay = response_state$sourceAssay,
        response_source_level = response_state$sourceLevel,
        lib_size = offset_state$libSize,
        block_spec = block_spec,
        genome_assay = offset_state$genomeAssay,
        offset_pseudocount = offset_state$offsetPseudocount
    )
}

#' Fit a Genome-Level Model
#'
#' `fitGenomeModel()` fits one negative-binomial model per genome using
#' `glmmTMB`. The response can come from a stored genome-level RNA assay or by
#' aggregating a gene-level assay to genomes on the fly.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`.
#'
#' This workflow is intended for the question of which genomes respond across
#' conditions or are associated with a continuous variable.
#'
#' @inheritParams makeGenomeModelData
#' @param keepFits Logical; if `TRUE`, store the backend `glmmTMB` model objects
#'   in the returned `MTTKFit`.
#' @param BPPARAM A `BiocParallelParam` instance controlling genome-wise model
#'   fitting.
#'
#' @return An `MTTKFit` with one row per genome.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitGenomeModel(x, variable = "condition")
#'     fit
#'     significantResults(fit)
#' }
#'
#' @export
fitGenomeModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = NULL,
    libraryOffset = TRUE,
    genomeOffset = NULL,
    libSize = NULL,
    sampleBlock = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL,
    keepFits = FALSE,
    BPPARAM = BiocParallel::SerialParam()
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use fitGenomeModel().",
            call. = FALSE
        )
    }

    if (!is.logical(keepFits) || length(keepFits) != 1L || is.na(keepFits)) {
        stop("'keepFits' must be TRUE or FALSE.", call. = FALSE)
    }

    model_spec <- .normalize_model_spec(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        referenceLevels = referenceLevels,
        allow_formula = TRUE,
        allow_random_slope = FALSE
    )
    .require_model_spec_term(model_spec)

    model_data <- makeGenomeModelData(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels
    )

    observations <- as.data.frame(model_data)
    model_state <- S4Vectors::metadata(model_data)$mttk_genome_model
    specification <- model_state$specification
    genome_summary <- model_state$genomeSummary
    response_assay <- model_state$responseAssay
    genome_assay <- model_state$genomeAssay
    variable_info <- list(
        type = model_state$variableType,
        referenceLevel = model_state$referenceLevel,
        contrastLevel = model_state$contrastLevel,
        effectLabel = model_state$effectLabel,
        modelSpec = model_spec
    )
    genome_ids <- rownames(genome_summary)

    fitted_rows <- BiocParallel::bplapply(genome_ids, function(genome_id) {
        data_genome <- observations[observations$genome_id == genome_id, , drop = FALSE]
        .fit_one_genome_model(
            data = data_genome,
            genome_id = genome_id,
            genome_summary = genome_summary[genome_id, , drop = FALSE],
            variable_info = variable_info,
            specification = specification,
            keep_fits = keepFits
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$result)
    )
    rownames(results) <- genome_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(one_fit) one_fit$model),
            genome_ids
        )
    } else {
        list()
    }

    formula <- .gene_model_formula(
        model_spec = variable_info$modelSpec,
        specification = specification,
        sample_block = !is.na(model_state$sampleBlock) && model_state$sampleBlock != ""
    )
    info <- list(
        backend = "glmmTMB",
        model = "genome_model",
        featureType = "genome",
        featureIdColumn = "genome_id",
        variable = model_state$variable,
        variableType = model_state$variableType,
        referenceLevel = model_state$referenceLevel,
        contrastLevel = model_state$contrastLevel,
        effectLabel = model_state$effectLabel,
        testedTermInput = model_state$testedTermInput,
        testedTerm = model_state$testedTerm,
        testedVariable = model_state$testedVariable,
        fixedEffectsFormula = model_state$fixedFormula,
        specification = specification,
        libraryOffset = model_state$libraryOffset,
        genomeOffset = model_state$genomeOffset,
        responseAssay = response_assay,
        sourceAssay = model_state$sourceAssay,
        sourceLevel = model_state$sourceLevel,
        libSizeSource = model_state$libSizeSource,
        genomeAssay = genome_assay,
        offsetPseudocount = model_state$offsetPseudocount,
        family = "nbinom2",
        formula = paste(deparse(formula), collapse = " "),
        randomEffects = .random_effect_description(
            sample_block = !is.na(model_state$sampleBlock) && model_state$sampleBlock != ""
        ),
        sampleBlock = model_state$sampleBlock,
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info = info,
        models = stored_models
    )
}
