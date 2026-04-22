.build_gene_model_observations <- function(
    x,
    model_spec,
    specification,
    assay_name,
    lib_size,
    block_spec,
    genome_assay,
    offset_pseudocount
) {
    counts <- SummarizedExperiment::assay(x, assay_name, withDimnames = TRUE)
    row_data <- SummarizedExperiment::rowData(x)
    gene_ids <- rownames(counts)
    sample_ids <- colnames(counts)

    gene_summary <- S4Vectors::DataFrame(
        gene_id = gene_ids,
        genome_id = as.character(row_data$genome_id),
        row.names = gene_ids
    )
    if ("gene_name" %in% names(row_data)) {
        gene_summary$gene_name <- as.character(row_data$gene_name)
    }

    obs <- expand.grid(
        gene_id = gene_ids,
        sample_id = sample_ids,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    obs$rna_count <- as.numeric(as.vector(counts))

    matched_genes <- match(obs$gene_id, gene_ids)
    obs$genome_id <- as.character(gene_summary$genome_id[matched_genes])
    if ("gene_name" %in% names(gene_summary)) {
        obs$gene_name <- as.character(gene_summary$gene_name[matched_genes])
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
        matched_genomes <- match(obs$genome_id, rownames(genome_mat))
        matched_samples <- match(obs$sample_id, colnames(genome_mat))

        if (anyNA(matched_genomes) || anyNA(matched_samples)) {
            stop(
                "Every gene/sample observation must be present in the selected genome assay.",
                call. = FALSE
            )
        }

        obs$genome_abundance <- as.numeric(genome_mat[cbind(matched_genomes, matched_samples)])
        obs$genome_abundance_offset <- obs$genome_abundance + offset_pseudocount
    }

    out <- S4Vectors::DataFrame(obs, check.names = FALSE)
    S4Vectors::metadata(out)$mttk_gene_model <- list(
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
        sourceAssay = assay_name,
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
        geneSummary = gene_summary
    )

    out
}

.empty_gene_model_row <- function(gene_id, gene_summary, variable_info) {
    columns <- list(
        gene_id = gene_id,
        genome_id = as.character(gene_summary$genome_id[[1L]])
    )
    if ("gene_name" %in% names(gene_summary)) {
        columns$gene_name <- as.character(gene_summary$gene_name[[1L]])
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
                row.names = gene_id
            )
        )
    )
}

.fit_one_gene_model <- function(data, gene_id, gene_summary, variable_info, model, keep_fits) {
    row <- .empty_gene_model_row(
        gene_id = gene_id,
        gene_summary = gene_summary,
        variable_info = variable_info
    )
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- "No observations were available for the gene."
        return(list(result = row, model = NULL))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- "All gene-level counts were zero."
        return(list(result = row, model = NULL))
    }

    warning_messages <- character()
    formula <- .gene_model_formula(
        model_spec = variable_info$modelSpec,
        specification = model,
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

#' Build a Gene-Level Model Data Table
#'
#' `makeGeneModelData()` materializes the long-form observation table used by
#' [fitGeneModel()]. Each row corresponds to one gene/sample observation, with
#' gene-level RNA counts, sample-level covariates, and an optional
#' parent-genome abundance offset aligned and ready for model fitting.
#'
#' This workflow is intended for gene-level differential expression or
#' association analysis, where the question is which individual genes change
#' across samples rather than which functions are associated across genomes.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)` for the
#'   simple one-variable interface. The column must be numeric or a factor with
#'   at least two levels. For factors with three or more levels supply `term` to
#'   select which contrast to test. Supply exactly one of `variable` or `formula`.
#' @param formula Optional one-sided or two-sided fixed-effect formula for the
#'   sample-level covariates, for example `~ condition + pH` or
#'   `rna_count ~ condition + pH`. Offsets are added internally by MTTK.
#'   Supply exactly one of `variable` or `formula`.
#' @param term Optional character string naming the model term to test. Required
#'   when `variable` refers to a factor with three or more levels (supply the
#'   contrast name, e.g. `"conditionB"`) or when `formula` defines more than one
#'   tested term.
#' @param assay Gene-level assay name used as the RNA response.
#' @param libraryOffset Logical; if `TRUE`, include `offset(log(lib_size))`.
#' @param libSize Library-size offset specification. Use `NULL` to compute
#'   `colSums(rnaGeneCounts(x))`, a single `colData(x)` column name, or a
#'   numeric vector with one value per sample.
#' @param sampleBlock Optional sample-level blocking variable from
#'   `colData(x)`. When supplied, MTTK adds a random intercept
#'   `(1 | sample_block)` to account for repeated measures or paired samples.
#' @param genomeOffset Logical; if `TRUE`, include
#'   `offset(log(genome_abundance + offsetPseudocount))`. If `NULL` (the
#'   default), MTTK uses genome-abundance normalization when `genomeAssay` is
#'   available in `genomeExperiment(x)`.
#' @param genomeAssay Genome-level assay used when `genomeOffset = TRUE`.
#' @param offsetPseudocount Non-negative pseudocount added to genome abundance
#'   before log-offset calculation.
#' @param referenceLevels Optional named list or named character vector setting
#'   the reference levels of factor-like sample covariates before model fitting.
#'
#' @return An `S4Vectors::DataFrame` with one row per gene/sample observation.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' gene_model_data <- makeGeneModelData(x, variable = "condition")
#' gene_model_data
#'
#' @export
makeGeneModelData <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
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

    assay_name <- .normalize_analysis_assays(x, assay)
    if (length(assay_name) != 1L) {
        stop("'assay' must be a single gene-level assay name.", call. = FALSE)
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

    offset_state <- .prepare_offset_inputs(
        x = x,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        libSize = libSize,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount
    )

    .build_gene_model_observations(
        x = x,
        model_spec = model_spec,
        specification = offset_state$specification,
        assay_name = assay_name,
        lib_size = offset_state$libSize,
        block_spec = block_spec,
        genome_assay = offset_state$genomeAssay,
        offset_pseudocount = offset_state$offsetPseudocount
    )
}

#' Fit a Gene-Level Model
#'
#' `fitGeneModel()` fits one negative-binomial model per gene using `glmmTMB`.
#' This workflow is intended for the gene-level question of which individual
#' genes change across conditions or are associated with a continuous variable.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`.
#'
#' In this workflow the parent genome is handled through the optional
#' genome-abundance offset rather than as a random effect, because each gene is
#' permanently assigned to one genome.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)` for the
#'   simple one-variable interface. The column must be numeric or a factor with
#'   at least two levels. For factors with three or more levels supply `term` to
#'   select which contrast to test. Supply exactly one of `variable` or `formula`.
#' @param formula Optional one-sided or two-sided fixed-effect formula for the
#'   sample-level covariates, for example `~ condition + pH` or
#'   `rna_count ~ condition + pH`. Offsets are added internally by MTTK.
#'   Supply exactly one of `variable` or `formula`.
#' @param term Optional character string naming the model term to test. Required
#'   when `variable` refers to a factor with three or more levels (supply the
#'   contrast name, e.g. `"conditionB"`) or when `formula` defines more than one
#'   tested term.
#' @param assay Gene-level assay name used as the RNA response.
#' @param libraryOffset Logical; if `TRUE`, include `offset(log(lib_size))`.
#' @param libSize Library-size offset specification. Use `NULL` to compute
#'   `colSums(rnaGeneCounts(x))`, a single `colData(x)` column name, or a
#'   numeric vector with one value per sample.
#' @param sampleBlock Optional sample-level blocking variable from
#'   `colData(x)`. When supplied, MTTK adds a random intercept
#'   `(1 | sample_block)` to account for repeated measures or paired samples.
#' @param genomeOffset Logical; if `TRUE`, include
#'   `offset(log(genome_abundance + offsetPseudocount))`. If `NULL` (the
#'   default), MTTK uses genome-abundance normalization when `genomeAssay` is
#'   available in `genomeExperiment(x)`.
#' @param genomeAssay Genome-level assay used when `genomeOffset = TRUE`.
#' @param offsetPseudocount Non-negative pseudocount added to genome abundance
#'   before log-offset calculation.
#' @param referenceLevels Optional named list or named character vector setting
#'   the reference levels of factor-like sample covariates before model fitting.
#' @param keepFits Logical; if `TRUE`, store the backend `glmmTMB` model objects
#'   in the returned `MTTKFit`.
#'
#' @return An `MTTKFit` with one row per gene.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitGeneModel(x, variable = "condition")
#'     fit
#'     significantResults(fit)
#' }
#'
#' @export
fitGeneModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
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
            "The 'glmmTMB' package must be installed to use fitGeneModel().",
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

    model_data <- makeGeneModelData(
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
    model_state <- S4Vectors::metadata(model_data)$mttk_gene_model
    specification <- model_state$specification
    gene_summary <- model_state$geneSummary
    assay_name <- model_state$sourceAssay
    genome_assay <- model_state$genomeAssay
    variable_info <- list(
        type = model_state$variableType,
        referenceLevel = model_state$referenceLevel,
        contrastLevel = model_state$contrastLevel,
        effectLabel = model_state$effectLabel,
        modelSpec = model_spec
    )
    gene_ids <- rownames(gene_summary)

    fitted_rows <- BiocParallel::bplapply(gene_ids, function(gene_id) {
        data_gene <- observations[observations$gene_id == gene_id, , drop = FALSE]
        .fit_one_gene_model(
            data = data_gene,
            gene_id = gene_id,
            gene_summary = gene_summary[gene_id, , drop = FALSE],
            variable_info = variable_info,
            model = specification,
            keep_fits = keepFits
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$result)
    )
    rownames(results) <- gene_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(one_fit) one_fit$model),
            gene_ids
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
        model = "gene_model",
        featureType = "gene",
        featureIdColumn = "gene_id",
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
        responseAssay = assay_name,
        sourceAssay = model_state$sourceAssay,
        sourceLevel = "gene",
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
