#' Find an mttk Workflow
#'
#' `findWorkflow()` is a static decision helper that maps a biological
#' question to the most relevant mttk workflow.
#'
#' The function is useful when you know the main unit of interest, such as
#' genes, genomes, KOs, modules, or pathways, but want help deciding which
#' modeling function best matches the question.
#'
#' The recommendations are intentionally question-focused:
#'
#' - `association`: an overall association or condition effect
#' - `heterogeneity`: different responses across genomes
#' - `group_difference`: a direct contrast between genome groups
#' - `group_coherence`: whether grouped genomes respond coherently
#' - `coherent_effects`: coherent KO-level effects within a module or pathway
#' - `total_activity`: total RNA assigned to a module or pathway
#' - `phylogenetic_signal`: whether fitted responses follow the tree globally
#' - `clade_shift`: whether a subtree shows a shifted response
#' - `phylogenetic_mean`: a phylogenetically corrected mean response
#'
#' @param level Primary unit of interest.
#' @param goal Biological question to answer.
#' @param variableType Type of predictor setup. Use `"categorical"` or
#'   `"continuous"` for the simple `variable =` interface, or `"formula"` when
#'   the model should adjust for additional covariates and test one focal term.
#' @param repeatedMeasures Logical; if `TRUE`, include a reminder that
#'   `sampleBlock` may be needed.
#' @param genomeOffset Logical or `NULL`. If `TRUE` or `FALSE`, include a
#'   reminder about explicitly setting `genomeOffset`. If `NULL`, no explicit
#'   offset recommendation is added.
#' @param phylogeny Logical; if `TRUE`, include a reminder about
#'   phylogeny-aware count models or tree-based follow-up where supported.
#'
#' @return An object of class `MTTKWorkflowRecommendation`, containing:
#'
#' - `question`: a short summary of the inferred biological question
#' - `recommendations`: an `S4Vectors::DataFrame` of primary and follow-up
#'   functions
#' - `suggestedArguments`: a character vector of arguments to consider
#' - `notes`: a character vector of workflow-specific notes
#'
#' @examples
#' findWorkflow(
#'     level = "ko",
#'     goal = "association",
#'     variableType = "categorical"
#' )
#'
#' findWorkflow(
#'     level = "module",
#'     goal = "coherent_effects",
#'     variableType = "continuous",
#'     phylogeny = TRUE
#' )
#'
#' @export
findWorkflow <- function(
    level = c("gene", "genome", "ko", "module", "pathway"),
    goal = c(
        "association",
        "heterogeneity",
        "group_difference",
        "group_coherence",
        "coherent_effects",
        "total_activity",
        "phylogenetic_signal",
        "clade_shift",
        "phylogenetic_mean"
    ),
    variableType = c("categorical", "continuous", "formula"),
    repeatedMeasures = FALSE,
    genomeOffset = NULL,
    phylogeny = FALSE
) {
    level <- match.arg(level)
    goal <- match.arg(goal)
    variableType <- match.arg(variableType)

    if (!is.logical(repeatedMeasures) || length(repeatedMeasures) != 1L || is.na(repeatedMeasures)) {
        stop("'repeatedMeasures' must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(genomeOffset) &&
        (!is.logical(genomeOffset) || length(genomeOffset) != 1L || is.na(genomeOffset))) {
        stop("'genomeOffset' must be TRUE, FALSE, or NULL.", call. = FALSE)
    }

    if (!is.logical(phylogeny) || length(phylogeny) != 1L || is.na(phylogeny)) {
        stop("'phylogeny' must be TRUE or FALSE.", call. = FALSE)
    }

    recommendation <- .recommend_mttk_workflow(level = level, goal = goal)
    suggested_arguments <- .workflow_argument_hints(
        variable_type = variableType,
        repeated_measures = repeatedMeasures,
        genome_offset = genomeOffset,
        phylogeny = phylogeny,
        level = level,
        goal = goal
    )

    notes <- unique(c(recommendation$notes, .workflow_extra_notes(
        level = level,
        goal = goal,
        phylogeny = phylogeny
    )))

    out <- list(
        question = recommendation$question,
        recommendations = recommendation$recommendations,
        suggestedArguments = suggested_arguments,
        notes = notes
    )
    class(out) <- c("MTTKWorkflowRecommendation", "list")
    out
}

.recommend_mttk_workflow <- function(level, goal) {
    make_steps <- function(functions, purposes) {
        S4Vectors::DataFrame(
            step = c("primary", rep("follow_up", length(functions) - 1L)),
            workflow_function = functions,
            purpose = purposes,
            row.names = NULL
        )
    }

    invalid <- function() {
        stop(
            "No mttk workflow is currently defined for level = '",
            level,
            "' and goal = '",
            goal,
            "'.",
            call. = FALSE
        )
    }

    if (identical(level, "gene")) {
        if (!identical(goal, "association")) {
            invalid()
        }

        return(list(
            question = "Which individual genes change or associate with the variable?",
            recommendations = make_steps(
                "fitGeneModel",
                "Fit one negative-binomial model per gene."
            ),
            notes = character()
        ))
    }

    if (identical(level, "genome")) {
        if (identical(goal, "association")) {
            return(list(
                question = "Which genomes show an overall transcriptional response?",
                recommendations = make_steps(
                    "fitGenomeModel",
                    "Fit one negative-binomial model per genome."
                ),
                notes = character()
            ))
        }

        if (identical(goal, "group_coherence")) {
            return(list(
                question = "Do genomes from the same group respond coherently?",
                recommendations = make_steps(
                    c("fitGenomeModel", "fitGenomeGroupMetaAnalysis"),
                    c(
                        "Estimate genome-level responses.",
                        "Summarize those responses within a genome group such as domain or clade."
                    )
                ),
                notes = "Set `group = \"...\"` to a grouping column in `genomeData(x)`."
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(list(
                question = "Do closely related genomes tend to respond similarly?",
                recommendations = make_steps(
                    c("fitGenomeModel", "fitGenomePhylogeneticSignal"),
                    c(
                        "Estimate genome-level responses.",
                        "Test for global phylogenetic signal in those responses."
                    )
                ),
                notes = "This workflow uses the genome tree as a follow-up analysis."
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(list(
                question = "Is there a subtree of genomes with a shifted response?",
                recommendations = make_steps(
                    c("fitGenomeModel", "scanGenomeClades"),
                    c(
                        "Estimate genome-level responses.",
                        "Scan the rooted tree for responding clades."
                    )
                ),
                notes = "A rooted `genomeTree(x)` is required for clade scans."
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(list(
                question = "What is the average genome response after phylogenetic correction?",
                recommendations = make_steps(
                    c("fitGenomeModel", "fitGenomePhylogeneticGLS"),
                    c(
                        "Estimate genome-level responses.",
                        "Estimate a phylogenetically corrected mean response."
                    )
                ),
                notes = "This workflow summarizes fitted genome responses with Brownian GLS."
            ))
        }

        invalid()
    }

    if (identical(level, "ko")) {
        if (identical(goal, "association")) {
            return(list(
                question = "Which KOs respond across genomes?",
                recommendations = make_steps(
                    "fitKOMixedModel",
                    "Fit one KO-level negative-binomial mixed model across genomes."
                ),
                notes = character()
            ))
        }

        if (identical(goal, "heterogeneity")) {
            return(list(
                question = "Does the same KO respond differently in different genomes?",
                recommendations = make_steps(
                    "fitKORandomSlopeModel",
                    "Fit one KO-level random-slope mixed model and estimate genome-specific KO responses."
                ),
                notes = "Use `koGenomeEffects()` to inspect the genome-specific conditional KO effects."
            ))
        }

        if (identical(goal, "group_difference")) {
            return(list(
                question = "Is the KO response directly different between genome groups?",
                recommendations = make_steps(
                    "fitKOGroupInteractionModel",
                    "Fit a direct KO-by-genome-group interaction model."
                ),
                notes = "This workflow requires a grouping column with exactly two usable groups."
            ))
        }

        if (identical(goal, "group_coherence")) {
            return(list(
                question = "Within genome groups, are KO responses coherent?",
                recommendations = make_steps(
                    c("fitKORandomSlopeModel", "fitKOGenomeGroupMetaAnalysis"),
                    c(
                        "Estimate genome-specific KO responses.",
                        "Summarize those KO responses within a genome group."
                    )
                ),
                notes = "Use `compareKOGenomeGroups()` afterward if you want a between-group contrast."
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(list(
                question = "Do genome-specific KO responses follow the phylogeny?",
                recommendations = make_steps(
                    c("fitKORandomSlopeModel", "fitKOPhylogeneticSignal"),
                    c(
                        "Estimate genome-specific KO responses.",
                        "Test whether those responses show global phylogenetic signal."
                    )
                ),
                notes = character()
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(list(
                question = "Is there a responding clade for a KO?",
                recommendations = make_steps(
                    c("fitKORandomSlopeModel", "scanKOClades"),
                    c(
                        "Estimate genome-specific KO responses.",
                        "Scan the rooted tree for clades with shifted KO responses."
                    )
                ),
                notes = "A rooted `genomeTree(x)` is required for clade scans."
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(list(
                question = "What is the average KO response after phylogenetic correction?",
                recommendations = make_steps(
                    c("fitKORandomSlopeModel", "fitKOPhylogeneticGLS"),
                    c(
                        "Estimate genome-specific KO responses.",
                        "Estimate a phylogenetically corrected KO-level mean response."
                    )
                ),
                notes = character()
            ))
        }

        invalid()
    }

    if (level %in% c("module", "pathway")) {
        feature_label <- level
        feature_title <- if (identical(level, "module")) "module" else "pathway"
        meta_fun <- if (identical(level, "module")) "fitModuleMetaAnalysis" else "fitPathwayMetaAnalysis"
        mixed_fun <- if (identical(level, "module")) "fitModuleMixedModel" else "fitPathwayMixedModel"
        random_fun <- if (identical(level, "module")) "fitModuleRandomSlopeModel" else "fitPathwayRandomSlopeModel"
        group_fun <- if (identical(level, "module")) "fitModuleGroupInteractionModel" else "fitPathwayGroupInteractionModel"
        signal_fun <- if (identical(level, "module")) "fitModulePhylogeneticSignal" else "fitPathwayPhylogeneticSignal"
        clade_fun <- if (identical(level, "module")) "scanModuleClades" else "scanPathwayClades"
        gls_fun <- if (identical(level, "module")) "fitModulePhylogeneticGLS" else "fitPathwayPhylogeneticGLS"

        if (identical(goal, "coherent_effects")) {
            return(list(
                question = paste0("Do the KOs assigned to this ", feature_title, " show a coherent direction of change?"),
                recommendations = make_steps(
                    c("fitKOMixedModel", meta_fun),
                    c(
                        "Estimate KO-level effects first.",
                        paste0("Summarize KO-level effects upward to the ", feature_title, ".")
                    )
                ),
                notes = c(
                    "For many-to-many KO annotations, `membershipMode = \"split\"` is usually the most conservative default.",
                    "This workflow summarizes KO effect sizes rather than refitting counts at the higher level."
                )
            ))
        }

        if (identical(goal, "total_activity")) {
            return(list(
                question = paste0("Does total RNA assigned to this ", feature_title, " change?"),
                recommendations = make_steps(
                    mixed_fun,
                    paste0("Re-aggregate counts to the ", feature_title, " level and refit the count model.")
                ),
                notes = c(
                    "`membershipMode = \"duplicate\"` targets total assigned activity.",
                    "`membershipMode = \"exclusive\"` restricts the analysis to uniquely assigned memberships."
                )
            ))
        }

        if (identical(goal, "heterogeneity")) {
            return(list(
                question = paste0("Does the same ", feature_title, " respond differently in different genomes?"),
                recommendations = make_steps(
                    random_fun,
                    paste0("Fit a random-slope model for genome-specific ", feature_title, " responses.")
                ),
                notes = character()
            ))
        }

        if (identical(goal, "group_difference")) {
            return(list(
                question = paste0("Is the ", feature_title, " response directly different between genome groups?"),
                recommendations = make_steps(
                    group_fun,
                    paste0("Fit a direct ", feature_title, "-by-genome-group interaction model.")
                ),
                notes = "This workflow requires a grouping column with exactly two usable groups."
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(list(
                question = paste0("Do genome-specific ", feature_title, " responses follow the phylogeny?"),
                recommendations = make_steps(
                    c(random_fun, signal_fun),
                    c(
                        paste0("Estimate genome-specific ", feature_title, " responses."),
                        paste0("Test whether those ", feature_title, " responses show global phylogenetic signal.")
                    )
                ),
                notes = character()
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(list(
                question = paste0("Is there a responding clade for a ", feature_title, "?"),
                recommendations = make_steps(
                    c(random_fun, clade_fun),
                    c(
                        paste0("Estimate genome-specific ", feature_title, " responses."),
                        paste0("Scan the rooted tree for clades with shifted ", feature_title, " responses.")
                    )
                ),
                notes = "A rooted `genomeTree(x)` is required for clade scans."
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(list(
                question = paste0("What is the average ", feature_title, " response after phylogenetic correction?"),
                recommendations = make_steps(
                    c(random_fun, gls_fun),
                    c(
                        paste0("Estimate genome-specific ", feature_title, " responses."),
                        paste0("Estimate a phylogenetically corrected mean ", feature_title, " response.")
                    )
                ),
                notes = character()
            ))
        }

        invalid()
    }

    invalid()
}

.workflow_argument_hints <- function(
    variable_type,
    repeated_measures,
    genome_offset,
    phylogeny,
    level,
    goal
) {
    hints <- if (identical(variable_type, "formula")) {
        c("formula = ~ ... + ...", "term = \"...\"")
    } else {
        "variable = \"...\""
    }

    if (isTRUE(repeated_measures)) {
        hints <- c(hints, "sampleBlock = \"...\"")
    }

    if (!is.null(genome_offset)) {
        hints <- c(hints, paste0("genomeOffset = ", genome_offset))
    }

    if (isTRUE(phylogeny) &&
        level %in% c("ko", "module", "pathway") &&
        goal %in% c("association", "heterogeneity", "group_difference", "total_activity")) {
        hints <- c(hints, "genomeCorrelation = \"brownian\"")
    }

    unique(hints)
}

.workflow_extra_notes <- function(level, goal, phylogeny) {
    notes <- character()

    if (isTRUE(phylogeny) &&
        level %in% c("ko", "module", "pathway") &&
        goal %in% c("association", "heterogeneity", "group_difference", "total_activity")) {
        notes <- c(
            notes,
            "When supported, `genomeCorrelation = \"brownian\"` uses the genome tree inside the count model."
        )
    }

    if (goal %in% c("phylogenetic_signal", "clade_shift", "phylogenetic_mean")) {
        notes <- c(notes, "These workflows require `genomeTree(x)` to be available.")
    }

    unique(notes)
}

#' @export
#' @noRd
print.MTTKWorkflowRecommendation <- function(x, ...) {
    cat("<MTTKWorkflowRecommendation>\n")
    cat("Question: ", x$question, "\n", sep = "")

    if (!is.null(x$recommendations) && nrow(x$recommendations) > 0L) {
        cat("\nRecommended functions:\n")
        recommendation_df <- as.data.frame(x$recommendations, stringsAsFactors = FALSE)
        for (i in seq_len(nrow(recommendation_df))) {
            cat(
                "- ",
                recommendation_df$step[[i]],
                ": ",
                recommendation_df$workflow_function[[i]],
                " - ",
                recommendation_df$purpose[[i]],
                "\n",
                sep = ""
            )
        }
    }

    if (!is.null(x$suggestedArguments) && length(x$suggestedArguments) > 0L) {
        cat("\nSuggested arguments:\n")
        for (one_arg in x$suggestedArguments) {
            cat("- ", one_arg, "\n", sep = "")
        }
    }

    if (!is.null(x$notes) && length(x$notes) > 0L) {
        cat("\nNotes:\n")
        for (one_note in x$notes) {
            cat("- ", one_note, "\n", sep = "")
        }
    }

    invisible(x)
}
