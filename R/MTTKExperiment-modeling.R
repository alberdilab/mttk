.normalize_ko_model_variable <- function(x, variable) {
    if (!is.character(variable) || length(variable) != 1L || is.na(variable) || variable == "") {
        stop("'variable' must be a single non-empty column name from colData(x).", call. = FALSE)
    }

    sample_data <- SummarizedExperiment::colData(x)
    if (!(variable %in% names(sample_data))) {
        stop("Unknown sample-level variable: ", variable, ".", call. = FALSE)
    }

    values <- sample_data[[variable]]
    sample_ids <- colnames(x)

    if (is.character(values)) {
        values <- factor(values, levels = unique(values))
    }

    if (is.factor(values)) {
        values <- droplevels(values)

        if (nlevels(values) != 2L) {
            stop(
                "'variable' must refer to a numeric column or a factor with exactly two levels.",
                call. = FALSE
            )
        }
    } else if (is.numeric(values) || is.integer(values)) {
        values <- as.numeric(values)

        if (length(unique(values)) < 2L) {
            stop(
                "'variable' must vary across samples before the model can be fit.",
                call. = FALSE
            )
        }
    } else {
        stop(
            "'variable' must refer to a numeric column or a factor/character column.",
            call. = FALSE
        )
    }

    if (anyNA(values)) {
        stop("'variable' must not contain missing values.", call. = FALSE)
    }

    stats::setNames(values, sample_ids)
}

.normalize_ko_model_lib_size <- function(x, libSize) {
    sample_ids <- colnames(x)

    if (is.null(libSize)) {
        values <- colSums(rnaGeneCounts(x))
        source <- "colSums(rna_gene_counts)"
    } else if (is.character(libSize) && length(libSize) == 1L && !is.na(libSize)) {
        sample_data <- SummarizedExperiment::colData(x)

        if (!(libSize %in% names(sample_data))) {
            stop("Unknown library-size column: ", libSize, ".", call. = FALSE)
        }

        values <- as.numeric(sample_data[[libSize]])
        source <- libSize
    } else if (is.numeric(libSize) && length(libSize) == ncol(x)) {
        values <- as.numeric(libSize)
        source <- "user_supplied"
    } else {
        stop(
            "'libSize' must be NULL, a single colData(x) column name, or a numeric vector with one value per sample.",
            call. = FALSE
        )
    }

    if (anyNA(values) || any(values <= 0)) {
        stop("'libSize' must contain strictly positive values.", call. = FALSE)
    }

    list(
        values = stats::setNames(values, sample_ids),
        source = source
    )
}

.gene_to_ko_link_df <- function(x) {
    feature_ids <- .normalize_feature_ids(x)
    link_list <- links(x)

    if ("gene_to_ko" %in% names(link_list)) {
        link_table <- as.data.frame(link_list[["gene_to_ko"]])

        if (!all(c("gene_id", "ko_id") %in% names(link_table))) {
            stop(
                "The 'gene_to_ko' link table must contain 'gene_id' and 'ko_id' columns.",
                call. = FALSE
            )
        }

        mapping <- unique(link_table[, c("gene_id", "ko_id"), drop = FALSE])
    } else {
        row_data <- SummarizedExperiment::rowData(x)

        if (!("ko_id" %in% names(row_data))) {
            stop(
                "A KO mapping is required in links(x)[['gene_to_ko']] or rowData(x)$ko_id.",
                call. = FALSE
            )
        }

        mapping <- data.frame(
            gene_id = feature_ids,
            ko_id = as.character(row_data$ko_id),
            stringsAsFactors = FALSE
        )
    }

    mapping <- mapping[!is.na(mapping$ko_id) & mapping$ko_id != "", , drop = FALSE]

    if (nrow(mapping) == 0L) {
        stop("No mapped KO identifiers were found for the modeled genes.", call. = FALSE)
    }

    if (!all(mapping$gene_id %in% feature_ids)) {
        stop(
            "The KO mapping contains gene identifiers that are not present in 'x'.",
            call. = FALSE
        )
    }

    mapping
}

.aggregate_rna_by_ko_genome <- function(x, assay_name) {
    gene_to_ko <- .gene_to_ko_link_df(x)
    gene_to_genome <- as.data.frame(.core_gene_to_genome_link(x))
    matched_genomes <- match(gene_to_ko$gene_id, gene_to_genome$gene_id)

    if (anyNA(matched_genomes)) {
        stop(
            "Every KO-mapped gene must also have a genome assignment before KO models can be fit.",
            call. = FALSE
        )
    }

    mapping <- data.frame(
        gene_id = gene_to_ko$gene_id,
        ko_id = gene_to_ko$ko_id,
        genome_id = as.character(gene_to_genome$genome_id[matched_genomes]),
        stringsAsFactors = FALSE
    )

    assay_mat <- SummarizedExperiment::assay(x, assay_name, withDimnames = TRUE)
    matched_genes <- match(mapping$gene_id, rownames(assay_mat))
    mapped_counts <- assay_mat[matched_genes, , drop = FALSE]
    pair_id <- paste(mapping$ko_id, mapping$genome_id, sep = "\r")
    pair_levels <- unique(pair_id)

    aggregated <- rowsum(
        mapped_counts,
        group = factor(pair_id, levels = pair_levels),
        reorder = FALSE
    )
    split_idx <- split(seq_len(nrow(mapping)), factor(pair_id, levels = pair_levels))
    pair_ids <- pair_levels

    pair_data <- S4Vectors::DataFrame(
        ko_id = vapply(split_idx, function(i) mapping$ko_id[i][1L], character(1)),
        genome_id = vapply(split_idx, function(i) mapping$genome_id[i][1L], character(1)),
        n_genes_pair = as.integer(vapply(
            split_idx,
            function(i) length(unique(mapping$gene_id[i])),
            integer(1)
        )),
        n_links_pair = as.integer(lengths(split_idx)),
        row.names = pair_ids
    )

    genome_data <- genomeData(x)
    matched_rows <- match(pair_data$genome_id, rownames(genome_data))
    if (all(!is.na(matched_rows))) {
        pair_data <- cbind(pair_data, genome_data[matched_rows, , drop = FALSE])
        rownames(pair_data) <- pair_ids
    }

    ko_ids <- unique(mapping$ko_id)
    ko_summary <- S4Vectors::DataFrame(
        ko_id = ko_ids,
        n_genes = as.integer(vapply(
            ko_ids,
            function(id) length(unique(mapping$gene_id[mapping$ko_id == id])),
            integer(1)
        )),
        n_links = as.integer(vapply(
            ko_ids,
            function(id) sum(mapping$ko_id == id),
            integer(1)
        )),
        n_genomes = as.integer(vapply(
            ko_ids,
            function(id) length(unique(mapping$genome_id[mapping$ko_id == id])),
            integer(1)
        )),
        row.names = ko_ids
    )

    list(
        mapping = mapping,
        counts = aggregated,
        pairData = pair_data,
        koSummary = ko_summary
    )
}

.build_ko_model_observations <- function(
    x,
    variable,
    model,
    assay_name,
    lib_size,
    genome_assay,
    offset_pseudocount
) {
    aggregated <- .aggregate_rna_by_ko_genome(x, assay_name = assay_name)
    pair_counts <- aggregated$counts
    pair_data <- aggregated$pairData
    sample_ids <- colnames(pair_counts)
    variable_values <- .normalize_ko_model_variable(x, variable)

    obs <- expand.grid(
        pair_id = rownames(pair_counts),
        sample_id = sample_ids,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    obs$rna_count <- as.numeric(as.vector(pair_counts))

    matched_pairs <- match(obs$pair_id, rownames(pair_data))
    obs$ko_id <- as.character(pair_data$ko_id[matched_pairs])
    obs$genome_id <- as.character(pair_data$genome_id[matched_pairs])
    obs$n_genes_pair <- as.integer(pair_data$n_genes_pair[matched_pairs])
    obs$n_links_pair <- as.integer(pair_data$n_links_pair[matched_pairs])

    obs[[variable]] <- variable_values[obs$sample_id]
    if (is.factor(variable_values)) {
        obs[[variable]] <- factor(
            obs[[variable]],
            levels = levels(variable_values)
        )
    }

    obs$lib_size <- as.numeric(lib_size$values[obs$sample_id])

    if (identical(model, "library_plus_genome_abundance")) {
        genome_mat <- SummarizedExperiment::assay(
            genomeExperiment(x),
            genome_assay,
            withDimnames = TRUE
        )
        matched_genomes <- match(as.character(pair_data$genome_id), rownames(genome_mat))

        if (anyNA(matched_genomes)) {
            stop(
                "Every KO/genome pair must be present in the selected genome assay.",
                call. = FALSE
            )
        }

        aligned_genome <- genome_mat[matched_genomes, sample_ids, drop = FALSE]
        obs$genome_abundance <- as.numeric(as.vector(aligned_genome))
        obs$genome_abundance_offset <- obs$genome_abundance + offset_pseudocount
    }

    list(
        observations = obs,
        koSummary = aggregated$koSummary
    )
}

.ko_model_formula <- function(variable, model) {
    rhs <- c(
        paste0("`", variable, "`"),
        "offset(log(lib_size))"
    )

    if (identical(model, "library_plus_genome_abundance")) {
        rhs <- c(rhs, "offset(log(genome_abundance_offset))")
    }

    rhs <- c(rhs, "(1 | genome_id)")
    stats::as.formula(paste("rna_count ~", paste(rhs, collapse = " + ")))
}

.empty_ko_model_row <- function(ko_id, ko_summary) {
    S4Vectors::DataFrame(
        ko_id = ko_id,
        tested_term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        q_value = NA_real_,
        intercept_estimate = NA_real_,
        intercept_std_error = NA_real_,
        n_genes = as.integer(ko_summary$n_genes[[1L]]),
        n_links = as.integer(ko_summary$n_links[[1L]]),
        n_genomes = as.integer(ko_summary$n_genomes[[1L]]),
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
        row.names = ko_id
    )
}

.fit_one_ko_mixed_model <- function(data, ko_id, ko_summary, variable, model, keep_fits) {
    row <- .empty_ko_model_row(ko_id, ko_summary)
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- "No observations were available for the KO."
        return(list(result = row, model = NULL))
    }

    if (length(unique(as.character(data$genome_id))) < 2L) {
        row$status <- "skipped"
        row$error_message <- "At least two genomes are required to estimate the random effect."
        return(list(result = row, model = NULL))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- "All KO-level counts were zero."
        return(list(result = row, model = NULL))
    }

    data$genome_id <- factor(as.character(data$genome_id))
    if (is.factor(data[[variable]])) {
        data[[variable]] <- droplevels(data[[variable]])
    }

    warning_messages <- character()
    formula <- .ko_model_formula(variable = variable, model = model)
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
    tested_terms <- setdiff(rownames(coefficient_table), "(Intercept)")

    if (length(tested_terms) != 1L) {
        row$status <- "error"
        row$error_message <- paste0(
            "Expected exactly one tested fixed-effect term, found ",
            length(tested_terms),
            "."
        )
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else {
            NA_character_
        }
        return(list(result = row, model = if (keep_fits) fit else NULL))
    }

    statistic_col <- intersect(colnames(coefficient_table), c("z value", "t value"))[1L]
    p_value_col <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)[1L]
    tested_row <- coefficient_table[tested_terms, , drop = FALSE]
    intercept_row <- coefficient_table["(Intercept)", , drop = FALSE]

    row$tested_term <- tested_terms
    row$estimate <- as.numeric(tested_row[, "Estimate"])
    row$std_error <- as.numeric(tested_row[, "Std. Error"])
    row$statistic <- as.numeric(tested_row[, statistic_col])
    row$p_value <- as.numeric(tested_row[, p_value_col])
    row$intercept_estimate <- as.numeric(intercept_row[, "Estimate"])
    row$intercept_std_error <- as.numeric(intercept_row[, "Std. Error"])
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

#' Fit a KO-Level Mixed Model
#'
#' `fitKOMixedModel()` fits one negative-binomial mixed model per KO using
#' `glmmTMB`. RNA counts are first aggregated to KO-within-genome observations,
#' so the fitted model respects the nesting of genes within genomes.
#'
#' The first implementation supports two model specifications:
#'
#' - `"library_only"` fits
#'   `rna_count ~ variable + offset(log(lib_size)) + (1 | genome_id)`
#' - `"library_plus_genome_abundance"` fits
#'   `rna_count ~ variable + offset(log(lib_size)) + offset(log(genome_abundance + offsetPseudocount)) + (1 | genome_id)`
#'
#' This workflow is intended for KO-level differential expression or association
#' analysis, where each KO can be observed in multiple genomes and the genome
#' term should be modeled explicitly as a random intercept.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)`. The
#'   column must be numeric or a factor with exactly two levels.
#' @param model Which of the two supported mixed-model specifications to fit.
#' @param assay Gene-level assay name used as the RNA response.
#' @param libSize Library-size offset specification. Use `NULL` to compute
#'   `colSums(rnaGeneCounts(x))`, a single `colData(x)` column name, or a
#'   numeric vector with one value per sample.
#' @param genomeAssay Genome-level assay used when
#'   `model = "library_plus_genome_abundance"`.
#' @param offsetPseudocount Non-negative pseudocount added to genome abundance
#'   before log-offset calculation.
#' @param keepFits Logical; if `TRUE`, store the backend `glmmTMB` model objects
#'   in the returned `MTTKFit`.
#'
#' @return An `MTTKFit` with one row per KO.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitKOMixedModel(x, variable = "condition")
#'     fit
#'     as.data.frame(fit)
#' }
#'
#' @export
fitKOMixedModel <- function(
    x,
    variable,
    model = c("library_only", "library_plus_genome_abundance"),
    assay = "rna_gene_counts",
    libSize = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    keepFits = FALSE
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use fitKOMixedModel().",
            call. = FALSE
        )
    }

    model <- match.arg(model)
    assay_name <- .normalize_analysis_assays(x, assay)
    if (length(assay_name) != 1L) {
        stop("'assay' must be a single gene-level assay name.", call. = FALSE)
    }
    lib_size <- .normalize_ko_model_lib_size(x, libSize = libSize)

    if (!is.logical(keepFits) || length(keepFits) != 1L || is.na(keepFits)) {
        stop("'keepFits' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.numeric(offsetPseudocount) ||
        length(offsetPseudocount) != 1L ||
        is.na(offsetPseudocount) ||
        offsetPseudocount < 0) {
        stop("'offsetPseudocount' must be a single non-negative numeric value.", call. = FALSE)
    }

    genome_assay <- if (identical(model, "library_plus_genome_abundance")) {
        .normalize_genome_analysis_assay(x, assay = genomeAssay)
    } else {
        NA_character_
    }

    model_data <- .build_ko_model_observations(
        x = x,
        variable = variable,
        model = model,
        assay_name = assay_name,
        lib_size = lib_size,
        genome_assay = genome_assay,
        offset_pseudocount = offsetPseudocount
    )

    observations <- model_data$observations
    ko_summary <- model_data$koSummary
    ko_ids <- rownames(ko_summary)

    fitted_rows <- lapply(ko_ids, function(ko_id) {
        data_ko <- observations[observations$ko_id == ko_id, , drop = FALSE]
        .fit_one_ko_mixed_model(
            data = data_ko,
            ko_id = ko_id,
            ko_summary = ko_summary[ko_id, , drop = FALSE],
            variable = variable,
            model = model,
            keep_fits = keepFits
        )
    })

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(x) x$result)
    )
    rownames(results) <- ko_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(x) x$model),
            ko_ids
        )
    } else {
        list()
    }

    formula <- .ko_model_formula(variable = variable, model = model)
    info <- list(
        backend = "glmmTMB",
        model = "ko_mixed_model",
        variable = variable,
        specification = model,
        responseAssay = assay_name,
        libSizeSource = lib_size$source,
        genomeAssay = genome_assay,
        offsetPseudocount = if (identical(model, "library_plus_genome_abundance")) {
            offsetPseudocount
        } else {
            NA_real_
        },
        family = "nbinom2",
        formula = paste(deparse(formula), collapse = " "),
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
