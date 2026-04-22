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
        n <- nlevels(values)

        if (n < 2L) {
            stop(
                "'variable' must refer to a factor with at least two levels.",
                call. = FALSE
            )
        }

        reference_level <- levels(values)[1L]

        if (n == 2L) {
            contrast_level <- levels(values)[2L]
            effect_label <- paste0(variable, ": ", contrast_level, " vs ", reference_level)
            variable_type <- "two_level_factor"
        } else {
            contrast_level <- NA_character_
            effect_label <- paste0(variable, ": contrasts vs ", reference_level)
            variable_type <- "multi_level_factor"
        }
    } else if (is.numeric(values) || is.integer(values)) {
        values <- as.numeric(values)

        if (length(unique(values)) < 2L) {
            stop(
                "'variable' must vary across samples before the model can be fit.",
                call. = FALSE
            )
        }

        reference_level <- NA_character_
        contrast_level <- NA_character_
        effect_label <- paste0(variable, " per unit increase")
        variable_type <- "numeric"
    } else {
        stop(
            "'variable' must refer to a numeric column or a factor/character column.",
            call. = FALSE
        )
    }

    if (anyNA(values)) {
        stop("'variable' must not contain missing values.", call. = FALSE)
    }

    list(
        name = variable,
        values = stats::setNames(values, sample_ids),
        type = variable_type,
        referenceLevel = reference_level,
        contrastLevel = contrast_level,
        effectLabel = effect_label
    )
}

.normalize_lib_size_values <- function(sample_data, sample_ids, default_values, default_label, libSize) {
    if (is.null(libSize)) {
        values <- as.numeric(default_values)
        source <- default_label
    } else if (is.character(libSize) && length(libSize) == 1L && !is.na(libSize)) {
        if (!(libSize %in% names(sample_data))) {
            stop("Unknown library-size column: ", libSize, ".", call. = FALSE)
        }

        values <- as.numeric(sample_data[[libSize]])
        source <- libSize
    } else if (is.numeric(libSize) && length(libSize) == length(sample_ids)) {
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

.normalize_ko_model_lib_size <- function(x, libSize) {
    .normalize_lib_size_values(
        sample_data = SummarizedExperiment::colData(x),
        sample_ids = colnames(x),
        default_values = colSums(rnaGeneCounts(x)),
        default_label = "colSums(rna_gene_counts)",
        libSize = libSize
    )
}

.normalize_model_variable <- .normalize_ko_model_variable
.normalize_model_lib_size <- .normalize_ko_model_lib_size

.normalize_reference_levels <- function(referenceLevels) {
    if (is.null(referenceLevels)) {
        return(list())
    }

    if (!is.list(referenceLevels)) {
        referenceLevels <- as.list(referenceLevels)
    }

    reference_names <- names(referenceLevels)
    if (is.null(reference_names) || anyNA(reference_names) || any(reference_names == "")) {
        stop(
            "'referenceLevels' must be NULL or a named list/vector of factor reference levels.",
            call. = FALSE
        )
    }

    normalized <- lapply(referenceLevels, function(x) {
        if (!is.character(x) || length(x) != 1L || is.na(x) || x == "") {
            stop(
                "Each entry of 'referenceLevels' must be a single non-empty character string.",
                call. = FALSE
            )
        }

        x
    })
    names(normalized) <- reference_names
    normalized
}

.coerce_model_column <- function(values, column_name, reference_level = NULL, require_binary_factor = FALSE) {
    if (is.character(values)) {
        values <- factor(values, levels = unique(values))
    } else if (is.logical(values)) {
        values <- factor(as.character(values), levels = unique(as.character(values)))
    } else if (is.factor(values)) {
        values <- droplevels(values)
    } else if (is.integer(values) || is.numeric(values)) {
        values <- as.numeric(values)
    } else {
        stop(
            "Modeling variable '",
            column_name,
            "' must be numeric, integer, logical, character, or factor.",
            call. = FALSE
        )
    }

    if (is.factor(values) && !is.null(reference_level)) {
        if (!(reference_level %in% levels(values))) {
            stop(
                "Reference level '",
                reference_level,
                "' is not present in sample variable '",
                column_name,
                "'.",
                call. = FALSE
            )
        }

        values <- stats::relevel(values, ref = reference_level)
        values <- droplevels(values)
    }

    if (anyNA(values)) {
        stop("Modeling variable '", column_name, "' must not contain missing values.", call. = FALSE)
    }

    if (is.factor(values)) {
        if (nlevels(values) < 2L) {
            stop(
                "Modeling factor '",
                column_name,
                "' must contain at least two levels.",
                call. = FALSE
            )
        }

    } else if (length(unique(values)) < 2L) {
        stop(
            "Modeling variable '",
            column_name,
            "' must vary across samples before the model can be fit.",
            call. = FALSE
        )
    }

    values
}

.model_sample_frame <- function(
    x,
    variables,
    referenceLevels = NULL,
    require_binary_variables = character()
) {
    sample_data <- as.data.frame(
        SummarizedExperiment::colData(x),
        stringsAsFactors = FALSE
    )
    sample_ids <- colnames(x)
    reference_levels <- .normalize_reference_levels(referenceLevels)

    variables <- unique(as.character(variables))
    if (length(variables) == 0L || anyNA(variables) || any(variables == "")) {
        stop("At least one valid sample-level modeling variable is required.", call. = FALSE)
    }

    missing_variables <- setdiff(variables, names(sample_data))
    if (length(missing_variables) > 0L) {
        stop(
            "Unknown sample-level variable(s): ",
            paste(missing_variables, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    frame <- sample_data[, variables, drop = FALSE]
    rownames(frame) <- sample_ids

    for (variable in variables) {
        frame[[variable]] <- .coerce_model_column(
            values = frame[[variable]],
            column_name = variable,
            reference_level = reference_levels[[variable]],
            require_binary_factor = variable %in% require_binary_variables
        )
    }

    frame
}

.normalize_sample_block_spec <- function(x, sampleBlock, fixed_variables = character()) {
    if (is.null(sampleBlock)) {
        return(list(
            name = NA_character_,
            values = NULL
        ))
    }

    if (!is.character(sampleBlock) ||
        length(sampleBlock) != 1L ||
        is.na(sampleBlock) ||
        sampleBlock == "") {
        stop(
            "'sampleBlock' must be NULL or a single non-empty column name from colData(x).",
            call. = FALSE
        )
    }

    if (sampleBlock %in% fixed_variables) {
        stop(
            "The sample-block variable '",
            sampleBlock,
            "' must not also be included among the fixed-effect variables.",
            call. = FALSE
        )
    }

    sample_data <- as.data.frame(
        SummarizedExperiment::colData(x),
        stringsAsFactors = FALSE
    )
    if (!(sampleBlock %in% names(sample_data))) {
        stop("Unknown sample-level block variable: ", sampleBlock, ".", call. = FALSE)
    }

    values <- sample_data[[sampleBlock]]
    if (is.factor(values)) {
        values <- droplevels(values)
    } else if (is.character(values) || is.logical(values) || is.numeric(values) || is.integer(values)) {
        values <- factor(as.character(values), levels = unique(as.character(values)))
    } else {
        stop(
            "The sample-block variable '",
            sampleBlock,
            "' must be character, factor, logical, integer, or numeric.",
            call. = FALSE
        )
    }

    if (anyNA(values) || any(levels(values) == "")) {
        stop(
            "The sample-block variable '",
            sampleBlock,
            "' must not contain missing or empty values.",
            call. = FALSE
        )
    }

    if (nlevels(values) < 2L) {
        stop(
            "The sample-block variable '",
            sampleBlock,
            "' must contain at least two levels.",
            call. = FALSE
        )
    }

    if (!any(table(values) > 1L)) {
        stop(
            "The sample-block variable '",
            sampleBlock,
            "' must repeat across at least two samples.",
            call. = FALSE
        )
    }

    list(
        name = sampleBlock,
        values = stats::setNames(values, colnames(x))
    )
}

.normalize_genome_correlation <- function(genomeCorrelation) {
    match.arg(genomeCorrelation, c("independent", "brownian"))
}

.resolve_feature_phylogeny <- function(object, genome_ids, genomeCorrelation, phylogeny = NULL) {
    genomeCorrelation <- .normalize_genome_correlation(genomeCorrelation)

    if (identical(genomeCorrelation, "independent")) {
        return(list(
            genomeCorrelation = genomeCorrelation,
            tree = NULL,
            treeSource = NA_character_,
            treeTipCount = NA_integer_
        ))
    }

    tree_source <- if (is.null(phylogeny)) "genomeTree(x)" else "phylogeny"
    tree <- .ensure_tree_branch_lengths(
        .require_genome_tree_object(object, tree = phylogeny)
    )

    genome_ids <- unique(as.character(genome_ids))
    genome_ids <- genome_ids[!is.na(genome_ids) & genome_ids != ""]
    missing_tips <- setdiff(genome_ids, tree$tip.label)
    if (length(missing_tips) > 0L) {
        stop(
            "The selected genome phylogeny is missing genome tip(s) required for this analysis: ",
            paste(missing_tips, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    list(
        genomeCorrelation = genomeCorrelation,
        tree = tree,
        treeSource = tree_source,
        treeTipCount = as.integer(length(tree$tip.label))
    )
}

.prepare_feature_genome_random_effect <- function(data, phylogeny_state) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)

    if (identical(phylogeny_state$genomeCorrelation, "brownian")) {
        feature_tree <- .prune_genome_tree(
            phylogeny_state$tree,
            unique(as.character(data$genome_id))
        )
        phylo_vcov <- ape::vcv(feature_tree)
        phylo_levels <- rownames(phylo_vcov)

        data$genome_id <- factor(as.character(data$genome_id), levels = phylo_levels)
        data$phylo_group <- factor("all_genomes", levels = "all_genomes")

        return(list(
            data = data,
            randomTerm = "propto(0 + genome_id | phylo_group, phylo_vcov)",
            formulaEnvironment = list(phylo_vcov = phylo_vcov)
        ))
    }

    data$genome_id <- factor(as.character(data$genome_id))
    list(
        data = data,
        randomTerm = "(1 | genome_id)",
        formulaEnvironment = NULL
    )
}

.set_formula_environment <- function(formula, values = NULL) {
    if (is.null(values) || length(values) == 0L) {
        return(formula)
    }

    formula_environment <- list2env(
        values,
        parent = environment(formula)
    )
    environment(formula) <- formula_environment
    formula
}

.fit_glmmtmb_nbinom2 <- function(formula, data, formula_environment = NULL) {
    formula <- .set_formula_environment(
        formula = formula,
        values = formula_environment
    )

    if (is.null(formula_environment) || length(formula_environment) == 0L) {
        return(
            glmmTMB::glmmTMB(
                formula = formula,
                data = data,
                family = glmmTMB::nbinom2(link = "log")
            )
        )
    }

    fit_environment <- list2env(
        c(
            list(
                formula = formula,
                data = data
            ),
            formula_environment
        ),
        parent = environment(formula)
    )

    eval(
        quote(
            glmmTMB::glmmTMB(
                formula = formula,
                data = data,
                family = glmmTMB::nbinom2(link = "log")
            )
        ),
        envir = fit_environment
    )
}

.normalize_fixed_effect_formula <- function(formula) {
    if (!inherits(formula, "formula")) {
        stop("'formula' must be a formula object.", call. = FALSE)
    }

    rhs <- if (length(formula) == 2L) {
        formula[[2L]]
    } else if (length(formula) == 3L) {
        formula[[3L]]
    } else {
        stop("'formula' must be a one-sided or two-sided formula.", call. = FALSE)
    }

    rhs_formula <- stats::as.formula(call("~", rhs))
    formula_text <- paste(deparse(rhs_formula), collapse = " ")

    if (grepl("\\|", formula_text)) {
        stop(
            "Random effects should not be included in 'formula'; MTTK adds the supported genome-level random effects internally.",
            call. = FALSE
        )
    }

    if (grepl("offset\\s*\\(", formula_text)) {
        stop(
            "Offsets should not be included in 'formula'; use 'libraryOffset' and 'genomeOffset' instead.",
            call. = FALSE
        )
    }

    rhs_formula
}

.available_model_terms <- function(fixed_formula, sample_data) {
    matrix <- stats::model.matrix(fixed_formula, data = sample_data)
    setdiff(colnames(matrix), "(Intercept)")
}

.term_matches_variable <- function(term, variable) {
    stripped_term <- gsub("`", "", as.character(term), fixed = TRUE)
    stripped_variable <- gsub("`", "", as.character(variable), fixed = TRUE)

    identical(stripped_term, stripped_variable) ||
        startsWith(stripped_term, stripped_variable)
}

.infer_tested_variable <- function(sample_data, tested_term) {
    if (is.null(tested_term) || is.na(tested_term)) {
        return(NA_character_)
    }

    variables <- names(sample_data)
    matches <- variables[vapply(
        variables,
        function(variable) .term_matches_variable(tested_term, variable),
        logical(1)
    )]

    if (length(matches) == 1L) {
        matches[[1L]]
    } else {
        NA_character_
    }
}

.describe_tested_term <- function(sample_data, tested_term, tested_variable = NA_character_) {
    tested_variable <- as.character(tested_variable)
    if (!is.na(tested_variable) &&
        tested_variable %in% names(sample_data)) {
        values <- sample_data[[tested_variable]]

        if (is.numeric(values)) {
            return(list(
                variable = tested_variable,
                type = "numeric",
                referenceLevel = NA_character_,
                contrastLevel = NA_character_,
                effectLabel = paste0(tested_variable, " per unit increase")
            ))
        }

        if (is.factor(values) && nlevels(values) == 2L) {
            return(list(
                variable = tested_variable,
                type = "two_level_factor",
                referenceLevel = levels(values)[1L],
                contrastLevel = levels(values)[2L],
                effectLabel = paste0(
                    tested_variable,
                    ": ",
                    levels(values)[2L],
                    " vs ",
                    levels(values)[1L]
                )
            ))
        }

        if (is.factor(values) && nlevels(values) >= 3L) {
            ref <- levels(values)[1L]
            stripped_term <- gsub("`", "", as.character(tested_term), fixed = TRUE)
            stripped_var  <- gsub("`", "", tested_variable, fixed = TRUE)
            contrast_level <- if (
                startsWith(stripped_term, stripped_var) &&
                nchar(stripped_term) > nchar(stripped_var)
            ) {
                substring(stripped_term, nchar(stripped_var) + 1L)
            } else {
                stripped_term
            }
            return(list(
                variable = tested_variable,
                type = "multi_level_factor",
                referenceLevel = ref,
                contrastLevel = contrast_level,
                effectLabel = paste0(tested_variable, ": ", contrast_level, " vs ", ref)
            ))
        }
    }

    list(
        variable = if (!is.na(tested_variable)) tested_variable else as.character(tested_term),
        type = "model_term",
        referenceLevel = NA_character_,
        contrastLevel = NA_character_,
        effectLabel = as.character(tested_term)
    )
}

.resolve_formula_random_slope <- function(sample_data, randomSlope, tested_variable, fixed_variables) {
    if (is.null(randomSlope)) {
        random_slope <- tested_variable
    } else {
        random_slope <- as.character(randomSlope)
    }

    if (length(random_slope) != 1L || is.na(random_slope) || random_slope == "") {
        return(NA_character_)
    }

    if (!(random_slope %in% fixed_variables)) {
        stop(
            "'randomSlope' must name a sample-level variable present in the fixed-effect formula.",
            call. = FALSE
        )
    }

    values <- sample_data[[random_slope]]
    if (!(is.numeric(values) || is.factor(values))) {
        stop(
            "The random-slope variable '",
            random_slope,
            "' must be numeric or factor-like.",
            call. = FALSE
        )
    }

    random_slope
}

.normalize_model_spec <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    referenceLevels = NULL,
    randomSlope = NULL,
    allow_formula = TRUE,
    allow_random_slope = FALSE
) {
    has_variable <- !is.null(variable)
    has_formula <- !is.null(formula)

    if (has_variable == has_formula) {
        stop("Provide exactly one of 'variable' or 'formula'.", call. = FALSE)
    }

    if (!is.null(term) &&
        (!is.character(term) || length(term) != 1L || is.na(term) || term == "")) {
        stop("'term' must be NULL or a single non-empty character string.", call. = FALSE)
    }

    if (has_variable) {
        variable <- as.character(variable)
        if (length(variable) != 1L || is.na(variable) || variable == "") {
            stop("'variable' must be a single non-empty column name from colData(x).", call. = FALSE)
        }

        sample_data <- .model_sample_frame(
            x = x,
            variables = variable,
            referenceLevels = referenceLevels
        )
        fixed_formula <- stats::as.formula(paste0("~ `", variable, "`"))
        candidate_terms <- .available_model_terms(fixed_formula, sample_data)
        if (is.null(term)) {
            if (length(candidate_terms) == 1L) {
                tested_term <- candidate_terms[[1L]]
            } else {
                stop(
                    "Variable '", variable, "' produces ", length(candidate_terms),
                    " contrasts (", paste(candidate_terms, collapse = ", "), "). ",
                    "Supply 'term' to select which contrast to test.",
                    call. = FALSE
                )
            }
        } else {
            tested_term <- .resolve_requested_model_term(
                term = term,
                candidates = candidate_terms,
                strict = TRUE,
                term_label = "term"
            )
        }
        term_info <- .describe_tested_term(
            sample_data = sample_data,
            tested_term = tested_term,
            tested_variable = variable
        )

        return(list(
            mode = "variable",
            sampleData = sample_data,
            fixedVariables = variable,
            fixedFormula = fixed_formula,
            fixedFormulaLabel = paste(deparse(fixed_formula), collapse = " "),
            availableTerms = candidate_terms,
            testedTermInput = if (is.null(term)) variable else term,
            testedTerm = tested_term,
            testedVariable = variable,
            variable = term_info$variable,
            variableType = term_info$type,
            referenceLevel = term_info$referenceLevel,
            contrastLevel = term_info$contrastLevel,
            effectLabel = term_info$effectLabel,
            randomSlopeVariable = if (allow_random_slope) variable else NA_character_
        ))
    }

    if (!allow_formula) {
        stop("This workflow does not support 'formula' yet.", call. = FALSE)
    }

    fixed_formula <- .normalize_fixed_effect_formula(formula)
    fixed_variables <- all.vars(fixed_formula)
    sample_data <- .model_sample_frame(
        x = x,
        variables = fixed_variables,
        referenceLevels = referenceLevels
    )
    candidate_terms <- .available_model_terms(fixed_formula, sample_data)
    tested_term <- if (length(candidate_terms) == 0L) {
        NA_character_
    } else if (is.null(term)) {
        if (length(candidate_terms) == 1L) candidate_terms[[1L]] else NA_character_
    } else {
        .resolve_requested_model_term(
            term = term,
            candidates = candidate_terms,
            strict = TRUE,
            term_label = "term"
        )
    }
    tested_variable <- .infer_tested_variable(sample_data, tested_term)
    term_info <- .describe_tested_term(
        sample_data = sample_data,
        tested_term = if (is.na(tested_term)) {
            if (!is.null(term)) term else NA_character_
        } else {
            tested_term
        },
        tested_variable = tested_variable
    )

    list(
        mode = "formula",
        sampleData = sample_data,
        fixedVariables = fixed_variables,
        fixedFormula = fixed_formula,
        fixedFormulaLabel = paste(deparse(fixed_formula), collapse = " "),
        availableTerms = candidate_terms,
        testedTermInput = if (!is.null(term)) term else NA_character_,
        testedTerm = tested_term,
        testedVariable = tested_variable,
        variable = term_info$variable,
        variableType = term_info$type,
        referenceLevel = term_info$referenceLevel,
        contrastLevel = term_info$contrastLevel,
        effectLabel = term_info$effectLabel,
        randomSlopeVariable = if (allow_random_slope) {
            .resolve_formula_random_slope(
                sample_data = sample_data,
                randomSlope = randomSlope,
                tested_variable = tested_variable,
                fixed_variables = fixed_variables
            )
        } else {
            NA_character_
        }
    )
}

.require_model_spec_term <- function(model_spec) {
    if ((is.null(model_spec$testedTerm) || is.na(model_spec$testedTerm)) &&
        length(model_spec$availableTerms) > 1L) {
        stop(
            "The fixed-effect formula produces multiple tested terms. Supply 'term' to select one of: ",
            paste(model_spec$availableTerms, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    invisible(model_spec)
}

.normalize_group_levels <- function(values, group_column, groupLevels = NULL) {
    values <- as.character(values)
    values <- values[!is.na(values) & values != ""]

    if (length(values) == 0L) {
        stop(
            "No usable genome group memberships were found in genomeData(x)$",
            group_column,
            ".",
            call. = FALSE
        )
    }

    if (is.null(groupLevels)) {
        levels <- unique(values)
        if (length(levels) != 2L) {
            stop(
                "Genome grouping column '",
                group_column,
                "' must contain exactly two usable groups for this workflow. ",
                "Supply 'groupLevels' to select two groups explicitly.",
                call. = FALSE
            )
        }

        return(levels)
    }

    if (!is.character(groupLevels) ||
        length(groupLevels) != 2L ||
        anyNA(groupLevels) ||
        any(groupLevels == "") ||
        length(unique(groupLevels)) != 2L) {
        stop(
            "'groupLevels' must be a character vector of length 2 giving the reference and contrast genome groups.",
            call. = FALSE
        )
    }

    missing_levels <- setdiff(groupLevels, unique(values))
    if (length(missing_levels) > 0L) {
        stop(
            "Unknown genome group level(s): ",
            paste(missing_levels, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    groupLevels
}

.normalize_feature_group_spec <- function(object, group, groupLevels = NULL) {
    mapping <- .genome_group_mapping_df(object, group_column = group)
    levels <- .normalize_group_levels(
        values = mapping$group_id,
        group_column = group,
        groupLevels = groupLevels
    )

    mapping <- mapping[mapping$group_id %in% levels, , drop = FALSE]
    mapping$genome_group <- factor(mapping$group_id, levels = levels)
    mapping <- mapping[!duplicated(mapping$genome_id), c("genome_id", "genome_group"), drop = FALSE]
    rownames(mapping) <- mapping$genome_id

    list(
        groupColumn = group,
        levels = levels,
        referenceLevel = levels[[1L]],
        contrastLevel = levels[[2L]],
        mapping = mapping
    )
}

.require_group_interaction_model_spec <- function(model_spec) {
    .require_model_spec_term(model_spec)

    if (is.null(model_spec$testedVariable) ||
        is.na(model_spec$testedVariable) ||
        model_spec$testedVariable == "") {
        stop(
            "Group-interaction models require a tested term that can be matched to a single sample-level variable.",
            call. = FALSE
        )
    }

    if (identical(model_spec$variableType, "model_term")) {
        stop(
            "Group-interaction models currently support numeric or two-level factor tested variables only.",
            call. = FALSE
        )
    }

    invisible(model_spec)
}

.build_feature_group_interaction_data <- function(model_data, object, group, groupLevels = NULL) {
    group_spec <- .normalize_feature_group_spec(
        object = object,
        group = group,
        groupLevels = groupLevels
    )

    obs <- as.data.frame(model_data, stringsAsFactors = FALSE)
    matched <- match(obs$genome_id, rownames(group_spec$mapping))
    keep <- !is.na(matched)
    obs <- obs[keep, , drop = FALSE]

    if (nrow(obs) == 0L) {
        stop(
            "No feature/genome observations could be aligned to the selected genome groups.",
            call. = FALSE
        )
    }

    obs$genome_group <- factor(
        as.character(group_spec$mapping$genome_group[matched[keep]]),
        levels = group_spec$levels
    )

    out <- S4Vectors::DataFrame(obs, check.names = FALSE)
    metadata_list <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    metadata_list$groupColumn <- group_spec$groupColumn
    metadata_list$groupReferenceLevel <- group_spec$referenceLevel
    metadata_list$groupContrastLevel <- group_spec$contrastLevel
    metadata_list$groupLevels <- group_spec$levels
    S4Vectors::metadata(out)$mttk_function_group_interaction <- metadata_list
    out
}

.feature_group_interaction_formula <- function(
    model_spec,
    specification,
    group_variable = "genome_group",
    genome_random_term = "(1 | genome_id)",
    sample_block = FALSE
) {
    tested_variable <- model_spec$testedVariable
    if (is.null(tested_variable) || is.na(tested_variable) || tested_variable == "") {
        stop(
            "Group-interaction models require a tested variable that can be identified from the fixed effects.",
            call. = FALSE
        )
    }

    rhs_terms <- c(
        paste(deparse(model_spec$fixedFormula[[2L]]), collapse = " "),
        paste0("`", group_variable, "`"),
        paste0("`", tested_variable, "`:`", group_variable, "`")
    )

    random_effects <- genome_random_term
    if (isTRUE(sample_block)) {
        random_effects <- c(random_effects, "(1 | sample_block)")
    }

    .compose_response_formula(
        fixed_formula = stats::as.formula(paste("~", paste(rhs_terms, collapse = " + "))),
        specification = specification,
        random_effects = random_effects
    )
}

.resolve_feature_group_interaction_term <- function(
    model_spec,
    coefficient_table,
    group_variable,
    group_contrast_level
) {
    tested_term <- .resolve_fitted_tested_term(model_spec, coefficient_table)
    candidates <- setdiff(rownames(coefficient_table), "(Intercept)")
    group_term <- .resolve_requested_model_term(
        term = paste0(group_variable, group_contrast_level),
        candidates = candidates,
        strict = TRUE,
        term_label = "group term"
    )

    strip_ticks <- function(x) {
        gsub("`", "", as.character(x), fixed = TRUE)
    }

    interaction_candidates <- candidates[grepl(":", candidates, fixed = TRUE)]
    tested_clean <- strip_ticks(tested_term)
    group_clean <- strip_ticks(group_term)
    matches <- interaction_candidates[vapply(interaction_candidates, function(candidate) {
        parts <- strsplit(strip_ticks(candidate), ":", fixed = TRUE)[[1L]]
        tested_clean %in% parts && group_clean %in% parts
    }, logical(1))]

    if (length(matches) != 1L) {
        stop(
            "Could not resolve a unique interaction coefficient for tested term '",
            tested_term,
            "' and genome group contrast '",
            group_contrast_level,
            "'.",
            call. = FALSE
        )
    }

    list(
        interactionTerm = matches[[1L]],
        baseTerm = tested_term,
        groupTerm = group_term
    )
}

.group_interaction_effect_label <- function(base_effect_label, group_spec) {
    paste0(
        base_effect_label,
        " difference in ",
        group_spec$contrastLevel,
        " vs ",
        group_spec$referenceLevel
    )
}

.empty_feature_group_interaction_row <- function(
    feature_id,
    feature_summary,
    variable_info,
    feature_id_column,
    group_spec
) {
    row <- .empty_feature_mixed_model_row(
        feature_id = feature_id,
        feature_summary = feature_summary,
        variable_info = variable_info,
        feature_id_column = feature_id_column
    )

    row$base_tested_term <- NA_character_
    row$group_term <- NA_character_
    row$group_column <- group_spec$groupColumn
    row$group_reference_level <- group_spec$referenceLevel
    row$group_contrast_level <- group_spec$contrastLevel
    row <- row[, c(
        feature_id_column,
        "base_tested_term",
        "tested_term",
        "group_term",
        "group_column",
        "group_reference_level",
        "group_contrast_level",
        setdiff(
            names(row),
            c(
                feature_id_column,
                "base_tested_term",
                "tested_term",
                "group_term",
                "group_column",
                "group_reference_level",
                "group_contrast_level"
            )
        )
    )]
    row
}

.fit_one_feature_group_interaction_model <- function(
    data,
    feature_id,
    feature_summary,
    model_spec,
    model,
    phylogeny_state,
    keep_fits,
    feature_id_column,
    feature_label,
    group_spec
) {
    variable_info <- list(
        type = paste0(model_spec$variableType, "_group_interaction"),
        referenceLevel = model_spec$referenceLevel,
        contrastLevel = model_spec$contrastLevel,
        effectLabel = .group_interaction_effect_label(
            base_effect_label = model_spec$effectLabel,
            group_spec = group_spec
        )
    )
    row <- .empty_feature_group_interaction_row(
        feature_id = feature_id,
        feature_summary = feature_summary,
        variable_info = variable_info,
        feature_id_column = feature_id_column,
        group_spec = group_spec
    )
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- paste0("No observations were available for the ", feature_label, ".")
        return(list(result = row, model = NULL))
    }

    if (length(unique(as.character(data$genome_id))) < 2L) {
        row$status <- "skipped"
        row$error_message <- "At least two genomes are required to estimate the random effect."
        return(list(result = row, model = NULL))
    }

    if (length(unique(as.character(data$genome_group))) < 2L) {
        row$status <- "skipped"
        row$error_message <- "At least two genome groups are required to estimate the interaction."
        return(list(result = row, model = NULL))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- paste0("All ", feature_label, "-level counts were zero.")
        return(list(result = row, model = NULL))
    }

    data$genome_group <- factor(as.character(data$genome_group), levels = group_spec$levels)
    random_effect_state <- .prepare_feature_genome_random_effect(
        data = data,
        phylogeny_state = phylogeny_state
    )
    data <- random_effect_state$data

    warning_messages <- character()
    formula <- .feature_group_interaction_formula(
        model_spec = model_spec,
        specification = model,
        genome_random_term = random_effect_state$randomTerm,
        sample_block = "sample_block" %in% names(data)
    )
    formula <- .set_formula_environment(
        formula = formula,
        values = random_effect_state$formulaEnvironment
    )
    fit <- tryCatch(
        withCallingHandlers(
            .fit_glmmtmb_nbinom2(
                formula = formula,
                data = data,
                formula_environment = random_effect_state$formulaEnvironment
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
    interaction_term <- tryCatch(
        .resolve_feature_group_interaction_term(
            model_spec = model_spec,
            coefficient_table = coefficient_table,
            group_variable = "genome_group",
            group_contrast_level = group_spec$contrastLevel
        ),
        error = identity
    )

    if (inherits(interaction_term, "error")) {
        row$status <- "error"
        row$error_message <- conditionMessage(interaction_term)
        row$warning_message <- if (length(warning_messages) > 0L) {
            paste(unique(warning_messages), collapse = " | ")
        } else {
            NA_character_
        }
        return(list(result = row, model = if (keep_fits) fit else NULL))
    }

    statistic_col <- intersect(colnames(coefficient_table), c("z value", "t value"))[1L]
    p_value_col <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)[1L]
    tested_row <- coefficient_table[interaction_term$interactionTerm, , drop = FALSE]
    intercept_row <- coefficient_table["(Intercept)", , drop = FALSE]
    row$base_tested_term <- interaction_term$baseTerm
    row$tested_term <- interaction_term$interactionTerm
    row$group_term <- interaction_term$groupTerm
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

.normalize_genome_model_lib_size <- function(x, response_counts, response_assay, libSize) {
    .normalize_lib_size_values(
        sample_data = SummarizedExperiment::colData(x),
        sample_ids = colnames(response_counts),
        default_values = colSums(response_counts),
        default_label = paste0("colSums(", response_assay, ")"),
        libSize = libSize
    )
}

.has_genome_analysis_assay <- function(x, assay) {
    if (!is.character(assay) || length(assay) != 1L || is.na(assay) || assay == "") {
        return(FALSE)
    }

    assay %in% names(genomeAssays(x, withDimnames = FALSE))
}

.offset_specification_name <- function(library_offset, genome_offset) {
    if (library_offset && genome_offset) {
        "library_plus_genome_abundance"
    } else if (library_offset) {
        "library_only"
    } else if (genome_offset) {
        "genome_only"
    } else {
        "no_offsets"
    }
}

.specification_uses_library_offset <- function(specification) {
    specification %in% c("library_only", "library_plus_genome_abundance")
}

.specification_uses_genome_offset <- function(specification) {
    specification %in% c("genome_only", "library_plus_genome_abundance")
}

.normalize_offset_flags <- function(x, libraryOffset, genomeOffset, genomeAssay) {
    if (!is.logical(libraryOffset) || length(libraryOffset) != 1L || is.na(libraryOffset)) {
        stop("'libraryOffset' must be TRUE or FALSE.", call. = FALSE)
    }

    if (is.null(genomeOffset)) {
        genomeOffset <- .has_genome_analysis_assay(x, genomeAssay)
    } else if (!is.logical(genomeOffset) || length(genomeOffset) != 1L || is.na(genomeOffset)) {
        stop("'genomeOffset' must be TRUE, FALSE, or NULL.", call. = FALSE)
    }

    list(
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        specification = .offset_specification_name(
            library_offset = libraryOffset,
            genome_offset = genomeOffset
        )
    )
}

.prepare_offset_inputs <- function(
    x,
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
        .normalize_ko_model_lib_size(x, libSize = libSize)
    } else {
        list(
            values = stats::setNames(rep(NA_real_, ncol(x)), colnames(x)),
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

.resolve_function_gene_mapping <- function(x, path) {
    path <- as.character(path)

    if (length(path) == 0L || anyNA(path) || any(path == "")) {
        stop("'path' must contain one or more named link tables.", call. = FALSE)
    }

    gene_ids <- rownames(x)
    if (is.null(gene_ids) || anyNA(gene_ids) || any(gene_ids == "")) {
        stop(
            "Row names must be present on 'x' before functional aggregation can be used.",
            call. = FALSE
        )
    }

    link_list <- links(x)
    available_links <- union(names(link_list), c("gene_to_genome", "gene_to_ko"))
    missing_links <- setdiff(path, available_links)

    if (length(missing_links) > 0L) {
        stop(
            "Unknown link table(s): ",
            paste(missing_links, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    mapping <- data.frame(
        gene_id = gene_ids,
        current_id = gene_ids,
        gene_index = seq_along(gene_ids),
        stringsAsFactors = FALSE
    )
    target_column <- "group_id"

    for (link_name in path) {
        link_table <- if (identical(link_name, "gene_to_ko")) {
            .gene_to_ko_link_df(x)
        } else if (identical(link_name, "gene_to_genome")) {
            .gene_to_genome_link_df(x)
        } else {
            as.data.frame(link_list[[link_name]])
        }

        if (ncol(link_table) < 2L) {
            stop(
                "Link table '",
                link_name,
                "' must contain at least two columns.",
                call. = FALSE
            )
        }

        target_column <- names(link_table)[2L]
        link_table <- unique(link_table[, 1:2, drop = FALSE])
        names(link_table) <- c("source_id", "target_id")

        mapping <- merge(
            mapping,
            link_table,
            by.x = "current_id",
            by.y = "source_id",
            all = FALSE,
            sort = FALSE
        )

        if (nrow(mapping) == 0L) {
            stop(
                "The link path did not map any genes: ",
                paste(path, collapse = " -> "),
                ".",
                call. = FALSE
            )
        }

        mapping <- mapping[
            order(mapping$gene_index),
            c("gene_id", "gene_index", "target_id"),
            drop = FALSE
        ]
        names(mapping)[3L] <- "current_id"
    }

    mapping <- unique(mapping[, c("gene_id", "gene_index", "current_id"), drop = FALSE])
    mapping <- mapping[order(mapping$gene_index), , drop = FALSE]
    names(mapping)[3L] <- "group_id"

    list(
        mapping = S4Vectors::DataFrame(
            gene_id = mapping$gene_id,
            group_id = mapping$group_id
        ),
        target_column = target_column,
        path = path
    )
}

.feature_genome_pair_id <- function(feature_id, genome_id) {
    if (length(feature_id) != length(genome_id)) {
        stop("'feature_id' and 'genome_id' must have the same length.", call. = FALSE)
    }

    paste(feature_id, genome_id, sep = "::")
}

.feature_summary_from_aggregation <- function(x, feature_id_column = NULL, feature_label = "feature") {
    state <- S4Vectors::metadata(x)$mttk_function_aggregation

    if (is.null(state) || is.null(state$featureSummary)) {
        stop(
            "The ",
            feature_label,
            "/genome aggregation does not contain stored feature summary metadata.",
            call. = FALSE
        )
    }

    if (!is.null(feature_id_column) &&
        !identical(state$featureIdColumn, feature_id_column)) {
        stop(
            "The aggregation stores '",
            state$featureIdColumn,
            "' summaries, not '",
            feature_id_column,
            "'.",
            call. = FALSE
        )
    }

    state$featureSummary
}

.feature_genome_observation_id <- function(feature_id_column) {
    paste0(sub("_id$", "", feature_id_column), "_genome_id")
}

.normalize_feature_membership_mode <- function(membership_mode) {
    match.arg(membership_mode, c("duplicate", "exclusive"))
}

.apply_feature_membership_mode <- function(mapping, membership_mode, feature_label) {
    membership_mode <- .normalize_feature_membership_mode(membership_mode)
    membership_count <- vapply(
        split(mapping$feature_id, mapping$gene_id),
        function(ids) length(unique(ids)),
        integer(1)
    )
    mapping$membership_count <- as.integer(membership_count[mapping$gene_id])

    if (identical(membership_mode, "exclusive")) {
        mapping <- mapping[mapping$membership_count == 1L, , drop = FALSE]
    }

    if (nrow(mapping) == 0L) {
        stop(
            "No usable ",
            feature_label,
            " mappings remained after applying membershipMode = '",
            membership_mode,
            "'.",
            call. = FALSE
        )
    }

    mapping
}

.aggregate_rna_by_feature_genome <- function(
    x,
    assay_name,
    path,
    feature_label,
    membership_mode = c("duplicate", "exclusive")
) {
    membership_mode <- .normalize_feature_membership_mode(membership_mode)
    resolved <- .resolve_function_gene_mapping(x, path)
    gene_to_genome <- as.data.frame(.core_gene_to_genome_link(x))
    matched_genomes <- match(resolved$mapping$gene_id, gene_to_genome$gene_id)

    if (anyNA(matched_genomes)) {
        stop(
            "Every ",
            feature_label,
            "-mapped gene must also have a genome assignment before mixed models can be fit.",
            call. = FALSE
        )
    }

    target_column <- resolved$target_column
    mapping <- data.frame(
        gene_id = as.character(resolved$mapping$gene_id),
        feature_id = as.character(resolved$mapping$group_id),
        genome_id = as.character(gene_to_genome$genome_id[matched_genomes]),
        stringsAsFactors = FALSE
    )
    input_link_count <- nrow(mapping)
    mapping <- .apply_feature_membership_mode(
        mapping = mapping,
        membership_mode = membership_mode,
        feature_label = feature_label
    )

    assay_mat <- SummarizedExperiment::assay(x, assay_name, withDimnames = TRUE)
    matched_genes <- match(mapping$gene_id, rownames(assay_mat))
    mapped_counts <- assay_mat[matched_genes, , drop = FALSE]
    pair_id <- .feature_genome_pair_id(mapping$feature_id, mapping$genome_id)
    pair_levels <- unique(pair_id)

    aggregated <- rowsum(
        mapped_counts,
        group = factor(pair_id, levels = pair_levels),
        reorder = FALSE
    )
    split_idx <- split(seq_len(nrow(mapping)), factor(pair_id, levels = pair_levels))
    pair_ids <- pair_levels

    pair_feature_ids <- vapply(
        split_idx,
        function(i) mapping$feature_id[i][1L],
        character(1)
    )
    pair_data <- S4Vectors::DataFrame(
        genome_id = vapply(split_idx, function(i) mapping$genome_id[i][1L], character(1)),
        n_genes_pair = as.integer(vapply(
            split_idx,
            function(i) length(unique(mapping$gene_id[i])),
            integer(1)
        )),
        n_links_pair = as.integer(lengths(split_idx)),
        row.names = pair_ids
    )
    pair_data[[target_column]] <- pair_feature_ids
    pair_data <- pair_data[, c(target_column, "genome_id", "n_genes_pair", "n_links_pair")]
    rownames(pair_data) <- pair_ids

    genome_data <- genomeData(x)
    matched_rows <- match(pair_data$genome_id, rownames(genome_data))
    if (nrow(genome_data) > 0L && all(!is.na(matched_rows))) {
        extra_columns <- setdiff(colnames(genome_data), colnames(pair_data))

        if (length(extra_columns) > 0L) {
            pair_data <- cbind(
                pair_data,
                genome_data[matched_rows, extra_columns, drop = FALSE]
            )
            rownames(pair_data) <- pair_ids
        }
    }

    feature_ids <- unique(mapping$feature_id)
    feature_summary <- S4Vectors::DataFrame(
        n_genes = as.integer(vapply(
            feature_ids,
            function(id) length(unique(mapping$gene_id[mapping$feature_id == id])),
            integer(1)
        )),
        n_links = as.integer(vapply(
            feature_ids,
            function(id) sum(mapping$feature_id == id),
            integer(1)
        )),
        n_genomes = as.integer(vapply(
            feature_ids,
            function(id) length(unique(mapping$genome_id[mapping$feature_id == id])),
            integer(1)
        )),
        row.names = feature_ids
    )
    feature_summary[[target_column]] <- feature_ids
    feature_summary <- feature_summary[, c(target_column, "n_genes", "n_links", "n_genomes")]
    rownames(feature_summary) <- feature_ids

    out <- SummarizedExperiment::SummarizedExperiment(
        assays = stats::setNames(list(aggregated), assay_name),
        rowData = pair_data,
        colData = SummarizedExperiment::colData(x)
    )
    S4Vectors::metadata(out)$mttk_function_aggregation <- list(
        sourceAssay = assay_name,
        featureSummary = feature_summary,
        featureIdColumn = target_column,
        featureLabel = feature_label,
        path = resolved$path,
        membershipMode = membership_mode,
        nInputLinks = as.integer(input_link_count),
        nRetainedLinks = as.integer(nrow(mapping))
    )
    out
}

.build_feature_mixed_model_observations <- function(
    aggregated,
    x,
    model_spec,
    specification,
    lib_size,
    block_spec,
    genome_assay,
    offset_pseudocount
) {
    aggregation_state <- S4Vectors::metadata(aggregated)$mttk_function_aggregation
    assay_name <- aggregation_state$sourceAssay
    feature_id_column <- aggregation_state$featureIdColumn
    observation_id_column <- .feature_genome_observation_id(feature_id_column)
    pair_counts <- SummarizedExperiment::assay(aggregated, assay_name, withDimnames = TRUE)
    pair_data <- SummarizedExperiment::rowData(aggregated)
    sample_ids <- colnames(pair_counts)

    obs <- expand.grid(
        feature_genome_id = rownames(pair_counts),
        sample_id = sample_ids,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    names(obs)[names(obs) == "feature_genome_id"] <- observation_id_column
    obs$rna_count <- as.numeric(as.vector(pair_counts))

    matched_pairs <- match(obs[[observation_id_column]], rownames(pair_data))
    obs[[feature_id_column]] <- as.character(pair_data[[feature_id_column]][matched_pairs])
    obs$genome_id <- as.character(pair_data$genome_id[matched_pairs])
    obs$n_genes_pair <- as.integer(pair_data$n_genes_pair[matched_pairs])
    obs$n_links_pair <- as.integer(pair_data$n_links_pair[matched_pairs])

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
        matched_genomes <- match(as.character(pair_data$genome_id), rownames(genome_mat))

        if (anyNA(matched_genomes)) {
            stop(
                "Every feature/genome pair must be present in the selected genome assay.",
                call. = FALSE
            )
        }

        aligned_genome <- genome_mat[matched_genomes, sample_ids, drop = FALSE]
        obs$genome_abundance <- as.numeric(as.vector(aligned_genome))
        obs$genome_abundance_offset <- obs$genome_abundance + offset_pseudocount
    }

    out <- S4Vectors::DataFrame(obs, check.names = FALSE)
    S4Vectors::metadata(out)$mttk_function_mixed_model <- list(
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
        randomSlopeVariable = model_spec$randomSlopeVariable,
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
        featureSummary = .feature_summary_from_aggregation(aggregated),
        featureIdColumn = feature_id_column,
        featureLabel = aggregation_state$featureLabel,
        path = aggregation_state$path,
        membershipMode = aggregation_state$membershipMode
    )

    out
}

.compose_response_formula <- function(fixed_formula, specification, random_effects = character()) {
    rhs <- paste(deparse(fixed_formula[[2L]]), collapse = " ")

    if (.specification_uses_library_offset(specification)) {
        rhs <- c(rhs, "offset(log(lib_size))")
    }

    if (.specification_uses_genome_offset(specification)) {
        rhs <- c(rhs, "offset(log(genome_abundance_offset))")
    }

    random_effects <- stats::na.omit(as.character(random_effects))
    random_effects <- random_effects[random_effects != ""]
    if (length(random_effects) > 0L) {
        rhs <- c(rhs, random_effects)
    }

    stats::as.formula(paste("rna_count ~", paste(rhs, collapse = " + ")))
}

.random_effect_description <- function(
    genome_random = FALSE,
    genome_slope = FALSE,
    genome_phylogenetic = FALSE,
    sample_block = FALSE
) {
    terms <- character()

    if (genome_slope) {
        if (isTRUE(genome_phylogenetic)) {
            terms <- c(terms, "genome_phylogenetic_random_intercept_and_slope")
        } else {
            terms <- c(terms, "genome_random_intercept_and_slope")
        }
    } else if (genome_random) {
        if (isTRUE(genome_phylogenetic)) {
            terms <- c(terms, "genome_phylogenetic_random_intercept")
        } else {
            terms <- c(terms, "genome_random_intercept")
        }
    }

    if (sample_block) {
        terms <- c(terms, "sample_block_random_intercept")
    }

    if (length(terms) == 0L) {
        NA_character_
    } else {
        paste(terms, collapse = " + ")
    }
}

.feature_model_formula <- function(
    model_spec,
    specification,
    genome_random_term = "(1 | genome_id)",
    sample_block = FALSE
) {
    random_effects <- genome_random_term
    if (isTRUE(sample_block)) {
        random_effects <- c(random_effects, "(1 | sample_block)")
    }

    .compose_response_formula(
        fixed_formula = model_spec$fixedFormula,
        specification = specification,
        random_effects = random_effects
    )
}

.feature_random_slope_formula <- function(model_spec, specification, sample_block = FALSE) {
    random_slope_variable <- model_spec$randomSlopeVariable

    if (is.na(random_slope_variable) || random_slope_variable == "") {
        stop(
            "A random-slope workflow requires a focal sample-level variable. Provide 'variable' or set 'randomSlope' when using 'formula'.",
            call. = FALSE
        )
    }

    random_effects <- paste0("(1 + `", random_slope_variable, "` | genome_id)")
    if (isTRUE(sample_block)) {
        random_effects <- c(random_effects, "(1 | sample_block)")
    }

    .compose_response_formula(
        fixed_formula = model_spec$fixedFormula,
        specification = specification,
        random_effects = random_effects
    )
}

.gene_model_formula <- function(model_spec, specification, sample_block = FALSE) {
    random_effects <- if (isTRUE(sample_block)) "(1 | sample_block)" else character()

    .compose_response_formula(
        fixed_formula = model_spec$fixedFormula,
        specification = specification,
        random_effects = random_effects
    )
}

.resolve_fitted_tested_term <- function(model_spec, coefficient_table) {
    tested_terms <- setdiff(rownames(coefficient_table), "(Intercept)")

    .resolve_requested_model_term(
        term = if (!is.null(model_spec$testedTerm) && !is.na(model_spec$testedTerm)) {
            model_spec$testedTerm
        } else {
            model_spec$testedTermInput
        },
        candidates = tested_terms,
        strict = TRUE,
        term_label = "term"
    )
}

.resolved_term_info <- function(model_spec, tested_term) {
    .describe_tested_term(
        sample_data = model_spec$sampleData,
        tested_term = tested_term,
        tested_variable = .infer_tested_variable(
            sample_data = model_spec$sampleData,
            tested_term = tested_term
        )
    )
}

.match_model_term_name <- function(term, candidates) {
    candidates <- as.character(candidates)
    if (term %in% candidates) {
        return(term)
    }

    strip_ticks <- function(x) {
        gsub("`", "", as.character(x), fixed = TRUE)
    }

    matched <- candidates[strip_ticks(candidates) == strip_ticks(term)]
    if (length(matched) == 1L) {
        return(matched[[1L]])
    }

    prefix_matches <- candidates[startsWith(strip_ticks(candidates), strip_ticks(term))]
    if (length(prefix_matches) == 1L) {
        return(prefix_matches[[1L]])
    }

    NA_character_
}

.resolve_requested_model_term <- function(term, candidates, strict = FALSE, term_label = "term") {
    candidates <- as.character(candidates)

    if (length(candidates) == 0L) {
        if (isTRUE(strict)) {
            stop("The fixed-effect formula does not define any tested terms.", call. = FALSE)
        }

        return(NA_character_)
    }

    if (is.null(term) || is.na(term) || term == "") {
        if (length(candidates) == 1L) {
            return(candidates[[1L]])
        }

        if (isTRUE(strict)) {
            stop(
                "The fixed-effect formula produces multiple tested terms. Supply '",
                term_label,
                "' to select one of: ",
                paste(candidates, collapse = ", "),
                ".",
                call. = FALSE
            )
        }

        return(NA_character_)
    }

    matched <- .match_model_term_name(term, candidates)
    if (is.na(matched) && isTRUE(strict)) {
        stop(
            "Could not match the requested '",
            term_label,
            "' to the available model terms. Available terms: ",
            paste(candidates, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    matched
}

.extract_glmmtmb_group_variance <- function(fit, group_name, tested_term) {
    vc <- tryCatch(
        glmmTMB::VarCorr(fit)$cond[[group_name]],
        error = function(e) NULL
    )

    if (is.null(vc)) {
        return(list(
            intercept_sd = NA_real_,
            slope_sd = NA_real_,
            intercept_slope_cor = NA_real_
        ))
    }

    stddev <- attr(vc, "stddev")
    if (is.null(stddev)) {
        return(list(
            intercept_sd = NA_real_,
            slope_sd = NA_real_,
            intercept_slope_cor = NA_real_
        ))
    }
    correlation <- attr(vc, "correlation")

    stddev_names <- names(stddev)
    if (is.null(stddev_names) && !is.null(dimnames(vc)[[1L]])) {
        stddev_names <- dimnames(vc)[[1L]]
    }
    names(stddev) <- stddev_names

    intercept_name <- .match_model_term_name("(Intercept)", stddev_names)
    slope_name <- .match_model_term_name(tested_term, stddev_names)

    intercept_sd <- if (!is.na(intercept_name) && intercept_name %in% names(stddev)) {
        as.numeric(stddev[[intercept_name]])
    } else {
        NA_real_
    }

    slope_sd <- if (!is.na(slope_name) && slope_name %in% names(stddev)) {
        as.numeric(stddev[[slope_name]])
    } else {
        NA_real_
    }

    intercept_slope_cor <- NA_real_
    if (!is.null(correlation)) {
        correlation_names <- rownames(correlation)
        intercept_cor_name <- .match_model_term_name("(Intercept)", correlation_names)
        slope_cor_name <- .match_model_term_name(tested_term, correlation_names)

        if (!is.na(intercept_cor_name) &&
            !is.na(slope_cor_name) &&
            intercept_cor_name %in% rownames(correlation) &&
            slope_cor_name %in% colnames(correlation)) {
            intercept_slope_cor <- as.numeric(correlation[intercept_cor_name, slope_cor_name])
        }
    }

    list(
        intercept_sd = intercept_sd,
        slope_sd = slope_sd,
        intercept_slope_cor = intercept_slope_cor
    )
}

.extract_glmmtmb_group_condvar <- function(fit, group_name) {
    ranef_fun <- get("ranef.glmmTMB", envir = asNamespace("glmmTMB"))
    group_effects <- tryCatch(
        ranef_fun(fit, condVar = TRUE)$cond[[group_name]],
        error = function(e) NULL
    )

    if (is.null(group_effects)) {
        return(NULL)
    }

    cond_var <- attr(group_effects, "condVar")
    list(
        effects = as.data.frame(group_effects, stringsAsFactors = FALSE),
        condVar = cond_var
    )
}

.extract_glmmtmb_fixed_effect_vcov <- function(fit) {
    vcov_mat <- tryCatch(
        as.matrix(stats::vcov(fit)$cond),
        error = function(e) NULL
    )

    if (is.null(vcov_mat) || !is.matrix(vcov_mat)) {
        return(NULL)
    }

    vcov_mat
}

.extract_glmmtmb_group_effects <- function(
    fit,
    feature_id,
    feature_id_column,
    group_name,
    tested_term,
    effect_label,
    fixed_intercept,
    fixed_effect
) {
    conditional_group <- tryCatch(
        stats::coef(fit)$cond[[group_name]],
        error = function(e) NULL
    )

    if (is.null(conditional_group)) {
        return(S4Vectors::DataFrame())
    }

    ranef_group <- .extract_glmmtmb_group_condvar(fit, group_name = group_name)
    fixed_vcov <- .extract_glmmtmb_fixed_effect_vcov(fit)

    conditional_group <- as.data.frame(conditional_group, stringsAsFactors = FALSE)
    conditional_group[[group_name]] <- rownames(conditional_group)

    conditional_term <- .match_model_term_name(tested_term, names(conditional_group))
    conditional_intercept <- .match_model_term_name("(Intercept)", names(conditional_group))

    group_ids <- as.character(conditional_group[[group_name]])
    conditional_intercept_estimate <- if (!is.na(conditional_intercept)) {
        as.numeric(conditional_group[[conditional_intercept]])
    } else {
        rep(NA_real_, nrow(conditional_group))
    }
    conditional_effect_estimate <- if (!is.na(conditional_term)) {
        as.numeric(conditional_group[[conditional_term]])
    } else {
        rep(NA_real_, nrow(conditional_group))
    }

    random_intercept_var <- rep(NA_real_, length(group_ids))
    random_effect_var <- rep(NA_real_, length(group_ids))
    if (!is.null(ranef_group) && !is.null(ranef_group$condVar)) {
        cond_var <- ranef_group$condVar
        random_term_names <- colnames(ranef_group$effects)

        intercept_random_name <- .match_model_term_name("(Intercept)", random_term_names)
        slope_random_name <- .match_model_term_name(tested_term, random_term_names)

        if (!is.na(intercept_random_name) && intercept_random_name %in% random_term_names) {
            intercept_index <- match(intercept_random_name, random_term_names)
            random_intercept_var <- vapply(
                seq_len(dim(cond_var)[3L]),
                function(i) as.numeric(cond_var[intercept_index, intercept_index, i]),
                numeric(1)
            )
        }

        if (!is.na(slope_random_name) && slope_random_name %in% random_term_names) {
            slope_index <- match(slope_random_name, random_term_names)
            random_effect_var <- vapply(
                seq_len(dim(cond_var)[3L]),
                function(i) as.numeric(cond_var[slope_index, slope_index, i]),
                numeric(1)
            )
        }
    }

    fixed_intercept_var <- NA_real_
    fixed_effect_var <- NA_real_
    if (!is.null(fixed_vcov)) {
        fixed_term_names <- rownames(fixed_vcov)
        fixed_intercept_name <- .match_model_term_name("(Intercept)", fixed_term_names)
        fixed_effect_name <- .match_model_term_name(tested_term, fixed_term_names)

        if (!is.na(fixed_intercept_name) && fixed_intercept_name %in% fixed_term_names) {
            fixed_intercept_var <- as.numeric(
                fixed_vcov[fixed_intercept_name, fixed_intercept_name]
            )
        }

        if (!is.na(fixed_effect_name) && fixed_effect_name %in% fixed_term_names) {
            fixed_effect_var <- as.numeric(
                fixed_vcov[fixed_effect_name, fixed_effect_name]
            )
        }
    }

    random_intercept_std_error <- sqrt(pmax(random_intercept_var, 0))
    random_effect_std_error <- sqrt(pmax(random_effect_var, 0))
    conditional_intercept_std_error <- sqrt(
        pmax(fixed_intercept_var + random_intercept_var, 0)
    )
    conditional_effect_std_error <- sqrt(
        pmax(fixed_effect_var + random_effect_var, 0)
    )

    out <- S4Vectors::DataFrame(
        tested_term = tested_term,
        effect_label = effect_label,
        conditional_intercept_estimate = conditional_intercept_estimate,
        conditional_intercept_std_error = conditional_intercept_std_error,
        conditional_effect_estimate = conditional_effect_estimate,
        conditional_effect_std_error = conditional_effect_std_error,
        random_intercept_std_error = random_intercept_std_error,
        random_effect_std_error = random_effect_std_error,
        random_intercept_deviation = conditional_intercept_estimate - fixed_intercept,
        random_effect_deviation = conditional_effect_estimate - fixed_effect,
        row.names = paste(feature_id, group_ids, sep = "::")
    )
    out[[feature_id_column]] <- feature_id
    out[[group_name]] <- group_ids
    out <- out[, c(
        feature_id_column,
        group_name,
        setdiff(names(out), c(feature_id_column, group_name))
    )]
    rownames(out) <- paste(feature_id, group_ids, sep = "::")
    out
}

.empty_feature_mixed_model_row <- function(feature_id, feature_summary, variable_info, feature_id_column) {
    columns <- list()
    columns[[feature_id_column]] <- feature_id

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
                n_genes = as.integer(feature_summary$n_genes[[1L]]),
                n_links = as.integer(feature_summary$n_links[[1L]]),
                n_genomes = as.integer(feature_summary$n_genomes[[1L]]),
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
                row.names = feature_id
            )
        )
    )
}

.fit_one_feature_mixed_model <- function(
    data,
    feature_id,
    feature_summary,
    model_spec,
    model,
    phylogeny_state,
    keep_fits,
    feature_id_column,
    feature_label
) {
    variable_info <- list(
        type = model_spec$variableType,
        referenceLevel = model_spec$referenceLevel,
        contrastLevel = model_spec$contrastLevel,
        effectLabel = model_spec$effectLabel
    )
    row <- .empty_feature_mixed_model_row(
        feature_id = feature_id,
        feature_summary = feature_summary,
        variable_info = variable_info,
        feature_id_column = feature_id_column
    )
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- paste0("No observations were available for the ", feature_label, ".")
        return(list(result = row, model = NULL))
    }

    if (length(unique(as.character(data$genome_id))) < 2L) {
        row$status <- "skipped"
        row$error_message <- "At least two genomes are required to estimate the random effect."
        return(list(result = row, model = NULL))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- paste0("All ", feature_label, "-level counts were zero.")
        return(list(result = row, model = NULL))
    }
    random_effect_state <- .prepare_feature_genome_random_effect(
        data = data,
        phylogeny_state = phylogeny_state
    )
    data <- random_effect_state$data

    warning_messages <- character()
    formula <- .feature_model_formula(
        model_spec = model_spec,
        specification = model,
        genome_random_term = random_effect_state$randomTerm,
        sample_block = "sample_block" %in% names(data)
    )
    formula <- .set_formula_environment(
        formula = formula,
        values = random_effect_state$formulaEnvironment
    )
    fit <- tryCatch(
        withCallingHandlers(
            .fit_glmmtmb_nbinom2(
                formula = formula,
                data = data,
                formula_environment = random_effect_state$formulaEnvironment
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
        .resolve_fitted_tested_term(model_spec, coefficient_table),
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
    intercept_row <- coefficient_table["(Intercept)", , drop = FALSE]
    term_info <- .resolved_term_info(model_spec, tested_term = tested_term)

    row$tested_term <- tested_term
    row$variable_type <- term_info$type
    row$reference_level <- term_info$referenceLevel
    row$contrast_level <- term_info$contrastLevel
    row$effect_label <- term_info$effectLabel
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

.empty_feature_random_slope_model_row <- function(
    feature_id,
    feature_summary,
    variable_info,
    feature_id_column
) {
    row <- .empty_feature_mixed_model_row(
        feature_id = feature_id,
        feature_summary = feature_summary,
        variable_info = variable_info,
        feature_id_column = feature_id_column
    )

    row$genome_intercept_sd <- NA_real_
    row$genome_slope_sd <- NA_real_
    row$genome_intercept_slope_cor <- NA_real_
    row
}

.fit_one_feature_random_slope_model <- function(
    data,
    feature_id,
    feature_summary,
    model_spec,
    model,
    keep_fits,
    feature_id_column,
    feature_label
) {
    variable_info <- list(
        type = model_spec$variableType,
        referenceLevel = model_spec$referenceLevel,
        contrastLevel = model_spec$contrastLevel,
        effectLabel = model_spec$effectLabel
    )
    row <- .empty_feature_random_slope_model_row(
        feature_id = feature_id,
        feature_summary = feature_summary,
        variable_info = variable_info,
        feature_id_column = feature_id_column
    )
    row$n_observations <- nrow(data)
    row$n_nonzero_observations <- sum(data$rna_count > 0)

    if (nrow(data) == 0L) {
        row$status <- "skipped"
        row$error_message <- paste0("No observations were available for the ", feature_label, ".")
        return(list(result = row, model = NULL, group_effects = S4Vectors::DataFrame()))
    }

    if (length(unique(as.character(data$genome_id))) < 2L) {
        row$status <- "skipped"
        row$error_message <- "At least two genomes are required to estimate random slopes."
        return(list(result = row, model = NULL, group_effects = S4Vectors::DataFrame()))
    }

    if (all(data$rna_count == 0)) {
        row$status <- "skipped"
        row$error_message <- paste0("All ", feature_label, "-level counts were zero.")
        return(list(result = row, model = NULL, group_effects = S4Vectors::DataFrame()))
    }

    data$genome_id <- factor(as.character(data$genome_id))

    warning_messages <- character()
    formula <- .feature_random_slope_formula(
        model_spec = model_spec,
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
        return(list(result = row, model = NULL, group_effects = S4Vectors::DataFrame()))
    }

    coefficient_table <- summary(fit)$coefficients$cond
    tested_term <- tryCatch(
        .resolve_fitted_tested_term(model_spec, coefficient_table),
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
        return(list(
            result = row,
            model = if (keep_fits) fit else NULL,
            group_effects = S4Vectors::DataFrame()
        ))
    }

    statistic_col <- intersect(colnames(coefficient_table), c("z value", "t value"))[1L]
    p_value_col <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)[1L]
    tested_row <- coefficient_table[tested_term, , drop = FALSE]
    intercept_row <- coefficient_table["(Intercept)", , drop = FALSE]
    variance_summary <- .extract_glmmtmb_group_variance(
        fit = fit,
        group_name = "genome_id",
        tested_term = tested_term
    )
    group_effects <- .extract_glmmtmb_group_effects(
        fit = fit,
        feature_id = feature_id,
        feature_id_column = feature_id_column,
        group_name = "genome_id",
        tested_term = tested_term,
        effect_label = .resolved_term_info(model_spec, tested_term)$effectLabel,
        fixed_intercept = as.numeric(intercept_row[, "Estimate"]),
        fixed_effect = as.numeric(tested_row[, "Estimate"])
    )
    term_info <- .resolved_term_info(model_spec, tested_term = tested_term)

    row$tested_term <- tested_term
    row$variable_type <- term_info$type
    row$reference_level <- term_info$referenceLevel
    row$contrast_level <- term_info$contrastLevel
    row$effect_label <- term_info$effectLabel
    row$estimate <- as.numeric(tested_row[, "Estimate"])
    row$std_error <- as.numeric(tested_row[, "Std. Error"])
    row$statistic <- as.numeric(tested_row[, statistic_col])
    row$p_value <- as.numeric(tested_row[, p_value_col])
    row$intercept_estimate <- as.numeric(intercept_row[, "Estimate"])
    row$intercept_std_error <- as.numeric(intercept_row[, "Std. Error"])
    row$genome_intercept_sd <- variance_summary$intercept_sd
    row$genome_slope_sd <- variance_summary$slope_sd
    row$genome_intercept_slope_cor <- variance_summary$intercept_slope_cor
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
        model = if (keep_fits) fit else NULL,
        group_effects = group_effects
    )
}

.fit_feature_random_slope_model <- function(
    x,
    variable,
    formula = NULL,
    term = NULL,
    assay,
    libraryOffset,
    genomeOffset,
    membershipMode,
    libSize,
    sampleBlock = NULL,
    genomeAssay,
    offsetPseudocount,
    referenceLevels = NULL,
    randomSlope = NULL,
    keepFits,
    path,
    feature_label,
    fit_model_name,
    BPPARAM = BiocParallel::SerialParam()
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use mixed-model workflows in MTTK.",
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
        randomSlope = randomSlope,
        allow_formula = TRUE,
        allow_random_slope = TRUE
    )
    .require_model_spec_term(model_spec)

    model_data <- .make_feature_mixed_model_data(
        x = x,
        model_spec = model_spec,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        path = path,
        feature_label = feature_label
    )

    observations <- as.data.frame(model_data)
    model_state <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    specification <- model_state$specification
    feature_summary <- model_state$featureSummary
    feature_id_column <- model_state$featureIdColumn
    assay_name <- model_state$sourceAssay
    genome_assay <- model_state$genomeAssay
    feature_ids <- rownames(feature_summary)

    fitted_rows <- BiocParallel::bplapply(feature_ids, function(feature_id) {
        data_feature <- observations[
            as.character(observations[[feature_id_column]]) == feature_id,
            ,
            drop = FALSE
        ]
        .fit_one_feature_random_slope_model(
            data = data_feature,
            feature_id = feature_id,
            feature_summary = feature_summary[feature_id, , drop = FALSE],
            model_spec = model_spec,
            model = specification,
            keep_fits = keepFits,
            feature_id_column = feature_id_column,
            feature_label = model_state$featureLabel
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$result)
    )
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(one_fit) one_fit$model),
            feature_ids
        )
    } else {
        list()
    }

    group_effects <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$group_effects)
    )
    if (is.null(group_effects)) {
        group_effects <- S4Vectors::DataFrame()
    }

    has_sample_block <- !is.na(model_state$sampleBlock) && model_state$sampleBlock != ""
    formula <- .feature_random_slope_formula(
        model_spec = model_spec,
        specification = specification,
        sample_block = has_sample_block
    )
    info <- list(
        backend = "glmmTMB",
        model = fit_model_name,
        featureType = model_state$featureLabel,
        featureIdColumn = feature_id_column,
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
        aggregationPath = model_state$path,
        membershipMode = model_state$membershipMode,
        libSizeSource = model_state$libSizeSource,
        genomeAssay = genome_assay,
        offsetPseudocount = model_state$offsetPseudocount,
        family = "nbinom2",
        formula = paste(deparse(formula), collapse = " "),
        randomEffects = .random_effect_description(
            genome_slope = TRUE,
            sample_block = has_sample_block
        ),
        groupEffectColumn = "genome_id",
        randomSlopeVariable = model_state$randomSlopeVariable,
        sampleBlock = model_state$sampleBlock,
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    out <- MTTKFit(
        results = results,
        info = info,
        models = stored_models
    )
    metadata_list <- S4Vectors::metadata(out)
    metadata_list$mttk_fit$groupEffects <- group_effects
    S4Vectors::metadata(out) <- metadata_list
    methods::validObject(out)
    out
}

.make_feature_mixed_model_data <- function(
    x,
    model_spec,
    assay,
    libraryOffset,
    genomeOffset,
    membershipMode,
    libSize,
    sampleBlock,
    genomeAssay,
    offsetPseudocount,
    path,
    feature_label
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    assay_name <- .normalize_analysis_assays(x, assay)
    if (length(assay_name) != 1L) {
        stop("'assay' must be a single gene-level assay name.", call. = FALSE)
    }

    offset_state <- .prepare_offset_inputs(
        x = x,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        libSize = libSize,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount
    )
    block_spec <- .normalize_sample_block_spec(
        x = x,
        sampleBlock = sampleBlock,
        fixed_variables = model_spec$fixedVariables
    )

    aggregated <- .aggregate_rna_by_feature_genome(
        x = x,
        assay_name = assay_name,
        path = path,
        feature_label = feature_label,
        membership_mode = membershipMode
    )
    .build_feature_mixed_model_observations(
        aggregated = aggregated,
        x = x,
        model_spec = model_spec,
        specification = offset_state$specification,
        lib_size = offset_state$libSize,
        block_spec = block_spec,
        genome_assay = offset_state$genomeAssay,
        offset_pseudocount = offset_state$offsetPseudocount
    )
}

.fit_feature_mixed_model <- function(
    x,
    variable,
    formula = NULL,
    term = NULL,
    assay,
    libraryOffset,
    genomeOffset,
    membershipMode,
    libSize,
    sampleBlock = NULL,
    genomeCorrelation = c("independent", "brownian"),
    phylogeny = NULL,
    genomeAssay,
    offsetPseudocount,
    referenceLevels = NULL,
    keepFits,
    path,
    feature_label,
    fit_model_name,
    BPPARAM = BiocParallel::SerialParam()
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use mixed-model workflows in MTTK.",
            call. = FALSE
        )
    }

    if (!is.logical(keepFits) || length(keepFits) != 1L || is.na(keepFits)) {
        stop("'keepFits' must be TRUE or FALSE.", call. = FALSE)
    }

    genome_correlation <- .normalize_genome_correlation(genomeCorrelation)

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

    model_data <- .make_feature_mixed_model_data(
        x = x,
        model_spec = model_spec,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        path = path,
        feature_label = feature_label
    )

    observations <- as.data.frame(model_data)
    model_state <- S4Vectors::metadata(model_data)$mttk_function_mixed_model
    specification <- model_state$specification
    feature_summary <- model_state$featureSummary
    feature_id_column <- model_state$featureIdColumn
    assay_name <- model_state$sourceAssay
    genome_assay <- model_state$genomeAssay
    feature_ids <- rownames(feature_summary)
    phylogeny_state <- .resolve_feature_phylogeny(
        object = x,
        genome_ids = observations$genome_id,
        genomeCorrelation = genome_correlation,
        phylogeny = phylogeny
    )

    fitted_rows <- BiocParallel::bplapply(feature_ids, function(feature_id) {
        data_feature <- observations[
            as.character(observations[[feature_id_column]]) == feature_id,
            ,
            drop = FALSE
        ]
        .fit_one_feature_mixed_model(
            data = data_feature,
            feature_id = feature_id,
            feature_summary = feature_summary[feature_id, , drop = FALSE],
            model_spec = model_spec,
            model = specification,
            phylogeny_state = phylogeny_state,
            keep_fits = keepFits,
            feature_id_column = feature_id_column,
            feature_label = model_state$featureLabel
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$result)
    )
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(one_fit) one_fit$model),
            feature_ids
        )
    } else {
        list()
    }

    has_sample_block <- !is.na(model_state$sampleBlock) && model_state$sampleBlock != ""
    genome_random_term <- if (identical(genome_correlation, "brownian")) {
        "propto(0 + genome_id | phylo_group, phylo_vcov)"
    } else {
        "(1 | genome_id)"
    }
    formula <- .feature_model_formula(
        model_spec = model_spec,
        specification = specification,
        genome_random_term = genome_random_term,
        sample_block = has_sample_block
    )
    info <- list(
        backend = "glmmTMB",
        model = fit_model_name,
        featureType = model_state$featureLabel,
        featureIdColumn = feature_id_column,
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
        aggregationPath = model_state$path,
        membershipMode = model_state$membershipMode,
        libSizeSource = model_state$libSizeSource,
        genomeAssay = genome_assay,
        offsetPseudocount = model_state$offsetPseudocount,
        genomeCorrelation = genome_correlation,
        treeSource = phylogeny_state$treeSource,
        treeTipCount = phylogeny_state$treeTipCount,
        family = "nbinom2",
        formula = paste(deparse(formula), collapse = " "),
        randomEffects = .random_effect_description(
            genome_random = TRUE,
            genome_phylogenetic = identical(genome_correlation, "brownian"),
            sample_block = has_sample_block
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

.fit_feature_group_interaction_model <- function(
    x,
    variable,
    formula = NULL,
    term = NULL,
    group,
    groupLevels = NULL,
    assay,
    libraryOffset,
    genomeOffset,
    membershipMode,
    libSize,
    sampleBlock = NULL,
    genomeCorrelation = c("independent", "brownian"),
    phylogeny = NULL,
    genomeAssay,
    offsetPseudocount,
    referenceLevels = NULL,
    keepFits,
    path,
    feature_label,
    fit_model_name,
    BPPARAM = BiocParallel::SerialParam()
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use mixed-model workflows in MTTK.",
            call. = FALSE
        )
    }

    if (!is.logical(keepFits) || length(keepFits) != 1L || is.na(keepFits)) {
        stop("'keepFits' must be TRUE or FALSE.", call. = FALSE)
    }

    genome_correlation <- .normalize_genome_correlation(genomeCorrelation)

    model_spec <- .normalize_model_spec(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        referenceLevels = referenceLevels,
        allow_formula = TRUE,
        allow_random_slope = FALSE
    )
    .require_group_interaction_model_spec(model_spec)

    model_data <- .make_feature_mixed_model_data(
        x = x,
        model_spec = model_spec,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        path = path,
        feature_label = feature_label
    )
    model_data <- .build_feature_group_interaction_data(
        model_data = model_data,
        object = x,
        group = group,
        groupLevels = groupLevels
    )

    observations <- as.data.frame(model_data)
    model_state <- S4Vectors::metadata(model_data)$mttk_function_group_interaction
    specification <- model_state$specification
    feature_summary <- model_state$featureSummary
    feature_id_column <- model_state$featureIdColumn
    assay_name <- model_state$sourceAssay
    genome_assay <- model_state$genomeAssay
    feature_ids <- rownames(feature_summary)
    phylogeny_state <- .resolve_feature_phylogeny(
        object = x,
        genome_ids = observations$genome_id,
        genomeCorrelation = genome_correlation,
        phylogeny = phylogeny
    )
    group_spec <- list(
        groupColumn = model_state$groupColumn,
        referenceLevel = model_state$groupReferenceLevel,
        contrastLevel = model_state$groupContrastLevel,
        levels = model_state$groupLevels
    )

    fitted_rows <- BiocParallel::bplapply(feature_ids, function(feature_id) {
        data_feature <- observations[
            as.character(observations[[feature_id_column]]) == feature_id,
            ,
            drop = FALSE
        ]
        .fit_one_feature_group_interaction_model(
            data = data_feature,
            feature_id = feature_id,
            feature_summary = feature_summary[feature_id, , drop = FALSE],
            model_spec = model_spec,
            model = specification,
            phylogeny_state = phylogeny_state,
            keep_fits = keepFits,
            feature_id_column = feature_id_column,
            feature_label = model_state$featureLabel,
            group_spec = group_spec
        )
    }, BPPARAM = BPPARAM)

    results <- do.call(
        rbind,
        lapply(fitted_rows, function(one_fit) one_fit$result)
    )
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    stored_models <- if (keepFits) {
        stats::setNames(
            lapply(fitted_rows, function(one_fit) one_fit$model),
            feature_ids
        )
    } else {
        list()
    }

    resolved_terms <- unique(stats::na.omit(as.character(results$tested_term)))

    has_sample_block <- !is.na(model_state$sampleBlock) && model_state$sampleBlock != ""
    genome_random_term <- if (identical(genome_correlation, "brownian")) {
        "propto(0 + genome_id | phylo_group, phylo_vcov)"
    } else {
        "(1 | genome_id)"
    }
    formula <- .feature_group_interaction_formula(
        model_spec = model_spec,
        specification = specification,
        genome_random_term = genome_random_term,
        sample_block = has_sample_block
    )
    info <- list(
        backend = "glmmTMB",
        model = fit_model_name,
        featureType = model_state$featureLabel,
        featureIdColumn = feature_id_column,
        variable = model_state$variable,
        variableType = paste0(model_state$variableType, "_group_interaction"),
        referenceLevel = model_state$referenceLevel,
        contrastLevel = model_state$contrastLevel,
        effectLabel = .group_interaction_effect_label(
            base_effect_label = model_state$effectLabel,
            group_spec = group_spec
        ),
        testedTermInput = model_state$testedTermInput,
        testedTerm = if (length(resolved_terms) == 1L) resolved_terms[[1L]] else NA_character_,
        testedBaseTerm = model_state$testedTerm,
        testedVariable = model_state$testedVariable,
        fixedEffectsFormula = model_state$fixedFormula,
        specification = specification,
        libraryOffset = model_state$libraryOffset,
        genomeOffset = model_state$genomeOffset,
        responseAssay = assay_name,
        aggregationPath = model_state$path,
        membershipMode = model_state$membershipMode,
        libSizeSource = model_state$libSizeSource,
        genomeAssay = genome_assay,
        offsetPseudocount = model_state$offsetPseudocount,
        genomeCorrelation = genome_correlation,
        treeSource = phylogeny_state$treeSource,
        treeTipCount = phylogeny_state$treeTipCount,
        family = "nbinom2",
        formula = paste(deparse(formula), collapse = " "),
        groupColumn = model_state$groupColumn,
        groupReferenceLevel = model_state$groupReferenceLevel,
        groupContrastLevel = model_state$groupContrastLevel,
        groupLevels = model_state$groupLevels,
        interactionVariable = model_state$testedVariable,
        randomEffects = .random_effect_description(
            genome_random = TRUE,
            genome_phylogenetic = identical(genome_correlation, "brownian"),
            sample_block = has_sample_block
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

.resolve_genome_model_response <- function(x, assay = NULL) {
    gene_assay_names <- names(SummarizedExperiment::assays(x, withDimnames = FALSE))
    genome_assay_names <- names(genomeAssays(x, withDimnames = FALSE))

    if (is.null(assay)) {
        assay <- if ("rna_genome_counts" %in% genome_assay_names) {
            "rna_genome_counts"
        } else if ("rna_gene_counts" %in% gene_assay_names) {
            "rna_gene_counts"
        } else {
            stop(
                "No default RNA assay was found. Provide a genome-level assay in genomeExperiment(x) or a gene-level assay in x.",
                call. = FALSE
            )
        }
    }

    assay <- as.character(assay)

    if (length(assay) != 1L || is.na(assay) || assay == "") {
        stop("'assay' must be NULL or a single non-empty assay name.", call. = FALSE)
    }

    if (assay %in% genome_assay_names) {
        genome_experiment <- genomeExperiment(x)
        response_counts <- SummarizedExperiment::assay(
            genome_experiment,
            assay,
            withDimnames = TRUE
        )
        genome_summary <- SummarizedExperiment::rowData(genome_experiment)
        if (!("genome_id" %in% names(genome_summary))) {
            genome_summary$genome_id <- rownames(genome_summary)
        }
        rownames(genome_summary) <- rownames(response_counts)

        return(list(
            counts = response_counts,
            genomeSummary = genome_summary,
            responseAssay = assay,
            sourceAssay = assay,
            sourceLevel = "genome_assay"
        ))
    }

    if (assay %in% gene_assay_names) {
        aggregated <- aggregateToGenome(x, assays = assay, fun = "sum")
        response_assay <- .to_genome_assay_name(assay)
        response_counts <- SummarizedExperiment::assay(
            aggregated,
            response_assay,
            withDimnames = TRUE
        )
        genome_summary <- SummarizedExperiment::rowData(aggregated)
        if (!("genome_id" %in% names(genome_summary))) {
            genome_summary$genome_id <- rownames(genome_summary)
        }
        rownames(genome_summary) <- rownames(response_counts)

        return(list(
            counts = response_counts,
            genomeSummary = genome_summary,
            responseAssay = response_assay,
            sourceAssay = assay,
            sourceLevel = "aggregated_from_gene_assay"
        ))
    }

    stop(
        "Unknown assay name: ",
        assay,
        ". Provide a gene-level assay in x or a genome-level assay in genomeExperiment(x).",
        call. = FALSE
    )
}

.ko_genome_pair_id <- function(ko_id, genome_id) {
    .feature_genome_pair_id(feature_id = ko_id, genome_id = genome_id)
}

.ko_summary_from_aggregation <- function(x) {
    .feature_summary_from_aggregation(
        x = x,
        feature_id_column = "ko_id",
        feature_label = "KO"
    )
}

.aggregate_rna_by_ko_genome <- function(x, assay_name) {
    out <- .aggregate_rna_by_feature_genome(
        x = x,
        assay_name = assay_name,
        path = "gene_to_ko",
        feature_label = "KO"
    )

    aggregation_state <- S4Vectors::metadata(out)$mttk_function_aggregation
    metadata_list <- S4Vectors::metadata(out)
    metadata_list$mttk_ko_aggregation <- list(
        sourceAssay = aggregation_state$sourceAssay,
        koSummary = aggregation_state$featureSummary
    )
    S4Vectors::metadata(out) <- metadata_list
    out
}

.build_ko_model_observations <- function(
    aggregated,
    x,
    model_spec,
    specification,
    lib_size,
    genome_assay,
    offset_pseudocount
) {
    .build_feature_mixed_model_observations(
        aggregated = aggregated,
        x = x,
        model_spec = model_spec,
        specification = specification,
        lib_size = lib_size,
        genome_assay = genome_assay,
        offset_pseudocount = offset_pseudocount
    )
}

.ko_model_formula <- function(model_spec, specification) {
    .feature_model_formula(model_spec = model_spec, specification = specification)
}

.empty_ko_model_row <- function(ko_id, ko_summary, variable_info) {
    .empty_feature_mixed_model_row(
        feature_id = ko_id,
        feature_summary = ko_summary,
        variable_info = variable_info,
        feature_id_column = "ko_id"
    )
}

.fit_one_ko_mixed_model <- function(
    data,
    ko_id,
    ko_summary,
    variable_info,
    model,
    keep_fits,
    phylogeny_state = list(genomeCorrelation = "independent", tree = NULL)
) {
    model_spec <- variable_info$modelSpec
    if (is.null(model_spec)) {
        stop(
            "Internal error: 'variable_info' must contain a 'modelSpec' entry.",
            call. = FALSE
        )
    }

    .fit_one_feature_mixed_model(
        data = data,
        feature_id = ko_id,
        feature_summary = ko_summary,
        model_spec = model_spec,
        model = model,
        phylogeny_state = phylogeny_state,
        keep_fits = keep_fits,
        feature_id_column = "ko_id",
        feature_label = "KO"
    )
}
