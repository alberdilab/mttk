#' Find an mttk Workflow
#'
#' `findWorkflow()` maps a biological question to the most relevant mttk
#' workflow.
#'
#' The function is designed as an interactive decision tree for regular users.
#' It steps through a small set of choices, such as the primary unit of
#' interest and the biological question, and then recommends the most relevant
#' mttk workflow.
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
#' @param ... Unused. `findWorkflow()` is interactive-only and does not accept
#'   direct workflow-selection arguments.
#'
#' @return An object of class `MTTKWorkflowRecommendation`, containing:
#'
#' - `question`: a short summary of the inferred biological question
#' - `recommendations`: an `S4Vectors::DataFrame` of primary and follow-up
#'   functions
#' - `suggestedArguments`: a character vector of arguments to consider
#' - `suggestedCode`: a character vector containing a ready-to-edit code
#'   template
#' - `notes`: a character vector of workflow-specific notes
#'
#' @examples
#' if (interactive()) {
#'     findWorkflow()
#' }
#'
#' @export
findWorkflow <- function(...) {
    dots <- list(...)

    if (length(dots) > 0L) {
        stop(
            "`findWorkflow()` is interactive-only and does not accept direct arguments. ",
            "Call `findWorkflow()` with no arguments and answer the prompts.",
            call. = FALSE
        )
    }

    if (!interactive()) {
        stop(
            "`findWorkflow()` is interactive-only. ",
            "Run it in an interactive R session and answer the prompts.",
            call. = FALSE
        )
    }

    selected <- .prompt_find_workflow(
        level = NULL,
        goal = NULL,
        variableType = NULL,
        repeatedMeasures = NULL,
        genomeOffset = NULL,
        phylogeny = NULL,
        promptLevel = TRUE,
        promptGoal = TRUE,
        promptVariableType = TRUE,
        promptRepeatedMeasures = TRUE,
        promptGenomeOffset = TRUE,
        promptPhylogeny = TRUE
    )

    .build_workflow_recommendation(
        level = selected$level,
        goal = selected$goal,
        variableType = selected$variableType,
        repeatedMeasures = selected$repeatedMeasures,
        genomeOffset = selected$genomeOffset,
        phylogeny = selected$phylogeny
    )
}

.build_workflow_recommendation <- function(
    level,
    goal,
    variableType,
    repeatedMeasures = FALSE,
    genomeOffset = NULL,
    phylogeny = FALSE
) {
    level <- .match_workflow_choice(level, .workflow_levels(), "level")
    goal <- .match_workflow_choice(goal, .workflow_goals(), "goal")
    variableType <- .match_workflow_choice(variableType, .workflow_variable_types(), "variableType")

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
    suggested_code <- .workflow_code_template(
        level = level,
        goal = goal,
        variableType = variableType,
        repeatedMeasures = repeatedMeasures,
        genomeOffset = genomeOffset,
        phylogeny = phylogeny
    )

    out <- list(
        question = recommendation$question,
        recommendations = recommendation$recommendations,
        suggestedArguments = suggested_arguments,
        suggestedCode = suggested_code,
        notes = notes
    )
    class(out) <- c("MTTKWorkflowRecommendation", "list")
    out
}

.workflow_predictor_args <- function(variableType) {
    switch(variableType,
        categorical = "variable = \"condition\"",
        continuous = "variable = \"oxygen\"",
        formula = c(
            "formula = ~ condition + salinity",
            "term = \"condition\""
        )
    )
}

.workflow_optional_args <- function(
    repeatedMeasures = FALSE,
    genomeOffset = NULL,
    phylogenyCountModel = FALSE,
    extraArgs = character()
) {
    args <- character()

    if (isTRUE(repeatedMeasures)) {
        args <- c(args, "sampleBlock = \"station_id\"")
    }

    if (!is.null(genomeOffset)) {
        args <- c(args, paste0("genomeOffset = ", genomeOffset))
    }

    if (isTRUE(phylogenyCountModel)) {
        args <- c(args, "genomeCorrelation = \"brownian\"")
    }

    c(args, extraArgs)
}

.workflow_call_lines <- function(assignTo, fun, positionalArgs, namedArgs = character()) {
    args <- c(positionalArgs, namedArgs)
    lines <- paste0("    ", args)

    if (length(lines) > 1L) {
        lines[-length(lines)] <- paste0(lines[-length(lines)], ",")
    }

    c(
        paste0(assignTo, " <- ", fun, "("),
        lines,
        ")"
    )
}

.workflow_code_template <- function(
    level,
    goal,
    variableType,
    repeatedMeasures,
    genomeOffset,
    phylogeny
) {
    predictor_args <- .workflow_predictor_args(variableType)
    count_phylogeny <- isTRUE(phylogeny) && .workflow_supports_phylo_count_model(level, goal)
    model_args <- c(
        predictor_args,
        .workflow_optional_args(
            repeatedMeasures = repeatedMeasures,
            genomeOffset = genomeOffset,
            phylogenyCountModel = count_phylogeny
        )
    )

    if (identical(level, "gene") && identical(goal, "association")) {
        return(.workflow_call_lines("gene_fit", "fitGeneModel", "x", model_args))
    }

    if (identical(level, "genome")) {
        if (identical(goal, "association")) {
            return(.workflow_call_lines("genome_fit", "fitGenomeModel", "x", model_args))
        }

        if (identical(goal, "group_coherence")) {
            return(c(
                .workflow_call_lines("genome_fit", "fitGenomeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "genome_group_fit",
                    "fitGenomeGroupMetaAnalysis",
                    c("genome_fit", "x"),
                    "group = \"domain\""
                )
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(c(
                .workflow_call_lines("genome_fit", "fitGenomeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "genome_signal_fit",
                    "fitGenomePhylogeneticSignal",
                    c("genome_fit", "x")
                )
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(c(
                .workflow_call_lines("genome_fit", "fitGenomeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "genome_clade_fit",
                    "scanGenomeClades",
                    c("genome_fit", "x")
                )
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(c(
                .workflow_call_lines("genome_fit", "fitGenomeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "genome_gls_fit",
                    "fitGenomePhylogeneticGLS",
                    c("genome_fit", "x")
                )
            ))
        }
    }

    if (identical(level, "ko")) {
        if (identical(goal, "association")) {
            return(.workflow_call_lines("ko_fit", "fitKOMixedModel", "x", model_args))
        }

        if (identical(goal, "heterogeneity")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKORandomSlopeModel", "x", model_args),
                "",
                "ko_effects <- koGenomeEffects(ko_fit)"
            ))
        }

        if (identical(goal, "group_difference")) {
            return(.workflow_call_lines(
                "ko_fit",
                "fitKOGroupInteractionModel",
                "x",
                c(model_args, "group = \"domain\"")
            ))
        }

        if (identical(goal, "group_coherence")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKORandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "ko_group_fit",
                    "fitKOGenomeGroupMetaAnalysis",
                    c("ko_fit", "x"),
                    "group = \"domain\""
                )
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKORandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "ko_signal_fit",
                    "fitKOPhylogeneticSignal",
                    c("ko_fit", "x")
                )
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKORandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "ko_clade_fit",
                    "scanKOClades",
                    c("ko_fit", "x")
                )
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKORandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "ko_gls_fit",
                    "fitKOPhylogeneticGLS",
                    c("ko_fit", "x")
                )
            ))
        }
    }

    if (identical(level, "module")) {
        if (identical(goal, "coherent_effects")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKOMixedModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "module_fit",
                    "fitModuleMetaAnalysis",
                    c("ko_fit", "x"),
                    "membershipMode = \"split\""
                )
            ))
        }

        if (identical(goal, "total_activity")) {
            return(.workflow_call_lines(
                "module_fit",
                "fitModuleMixedModel",
                "x",
                c(model_args, "membershipMode = \"duplicate\"")
            ))
        }

        if (identical(goal, "heterogeneity")) {
            return(c(
                .workflow_call_lines("module_fit", "fitModuleRandomSlopeModel", "x", model_args),
                "",
                "module_effects <- moduleGenomeEffects(module_fit)"
            ))
        }

        if (identical(goal, "group_difference")) {
            return(.workflow_call_lines(
                "module_fit",
                "fitModuleGroupInteractionModel",
                "x",
                c(model_args, "group = \"domain\"")
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(c(
                .workflow_call_lines("module_fit", "fitModuleRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "module_signal_fit",
                    "fitModulePhylogeneticSignal",
                    c("module_fit", "x")
                )
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(c(
                .workflow_call_lines("module_fit", "fitModuleRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "module_clade_fit",
                    "scanModuleClades",
                    c("module_fit", "x")
                )
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(c(
                .workflow_call_lines("module_fit", "fitModuleRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "module_gls_fit",
                    "fitModulePhylogeneticGLS",
                    c("module_fit", "x")
                )
            ))
        }
    }

    if (identical(level, "pathway")) {
        if (identical(goal, "coherent_effects")) {
            return(c(
                .workflow_call_lines("ko_fit", "fitKOMixedModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "pathway_fit",
                    "fitPathwayMetaAnalysis",
                    c("ko_fit", "x"),
                    "membershipMode = \"split\""
                )
            ))
        }

        if (identical(goal, "total_activity")) {
            return(.workflow_call_lines(
                "pathway_fit",
                "fitPathwayMixedModel",
                "x",
                c(model_args, "membershipMode = \"duplicate\"")
            ))
        }

        if (identical(goal, "heterogeneity")) {
            return(c(
                .workflow_call_lines("pathway_fit", "fitPathwayRandomSlopeModel", "x", model_args),
                "",
                "pathway_effects <- pathwayGenomeEffects(pathway_fit)"
            ))
        }

        if (identical(goal, "group_difference")) {
            return(.workflow_call_lines(
                "pathway_fit",
                "fitPathwayGroupInteractionModel",
                "x",
                c(model_args, "group = \"domain\"")
            ))
        }

        if (identical(goal, "phylogenetic_signal")) {
            return(c(
                .workflow_call_lines("pathway_fit", "fitPathwayRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "pathway_signal_fit",
                    "fitPathwayPhylogeneticSignal",
                    c("pathway_fit", "x")
                )
            ))
        }

        if (identical(goal, "clade_shift")) {
            return(c(
                .workflow_call_lines("pathway_fit", "fitPathwayRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "pathway_clade_fit",
                    "scanPathwayClades",
                    c("pathway_fit", "x")
                )
            ))
        }

        if (identical(goal, "phylogenetic_mean")) {
            return(c(
                .workflow_call_lines("pathway_fit", "fitPathwayRandomSlopeModel", "x", model_args),
                "",
                .workflow_call_lines(
                    "pathway_gls_fit",
                    "fitPathwayPhylogeneticGLS",
                    c("pathway_fit", "x")
                )
            ))
        }
    }

    character()
}

.workflow_levels <- function() {
    c("gene", "genome", "ko", "module", "pathway")
}

.workflow_goals <- function() {
    c(
        "association",
        "heterogeneity",
        "group_difference",
        "group_coherence",
        "coherent_effects",
        "total_activity",
        "phylogenetic_signal",
        "clade_shift",
        "phylogenetic_mean"
    )
}

.workflow_variable_types <- function() {
    c("categorical", "continuous", "formula")
}

.match_workflow_choice <- function(value, choices, arg) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || value == "") {
        stop("'", arg, "' must be a single non-empty character value.", call. = FALSE)
    }

    if (!(value %in% choices)) {
        stop(
            "'",
            arg,
            "' must be one of: ",
            paste(shQuote(choices), collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    value
}

.workflow_level_labels <- function() {
    c(
        gene = paste(
            "Gene: individual genes as the unit of inference",
            "  Example: which transporters or stress-response genes change between treatments?",
            sep = "\n"
        ),
        genome = paste(
            "Genome: whole-genome responses",
            "  Example: which MAGs shift their overall transcription across an oxygen gradient?",
            sep = "\n"
        ),
        ko = paste(
            "KO: shared functions across genomes",
            "  Example: which KOs respond across the community, regardless of which genome encodes them?",
            sep = "\n"
        ),
        module = paste(
            "Module: higher functional groups built from KOs",
            "  Example: do methanogenesis modules show a coordinated response?",
            sep = "\n"
        ),
        pathway = paste(
            "Pathway: broader functional pathways built from KOs",
            "  Example: does the total RNA assigned to sulfur metabolism pathways change?",
            sep = "\n"
        )
    )
}

.workflow_goal_labels <- function(level) {
    labels <- c(
        association = paste(
            "Overall association or condition effect",
            "  Example: which genes, genomes, or KOs increase in the treatment group?",
            sep = "\n"
        ),
        heterogeneity = paste(
            "Different responses across genomes",
            "  Example: does the same KO increase in some genomes but decrease in others?",
            sep = "\n"
        ),
        group_difference = paste(
            "Direct difference between genome groups",
            "  Example: does a KO respond differently in Archaea versus Bacteria?",
            sep = "\n"
        ),
        group_coherence = paste(
            "Coherent responses within genome groups",
            "  Example: do genomes from the same clade show a consistent overall response?",
            sep = "\n"
        ),
        coherent_effects = paste(
            "Coherent KO-level effects within the higher functional category",
            "  Example: do the KOs assigned to a pathway all point in the same direction of change?",
            sep = "\n"
        ),
        total_activity = paste(
            "Total RNA assigned to the higher functional category",
            "  Example: does the total activity of a module increase, regardless of KO-level agreement?",
            sep = "\n"
        ),
        phylogenetic_signal = paste(
            "Global phylogenetic signal in fitted responses",
            "  Example: do closely related genomes tend to respond similarly?",
            sep = "\n"
        ),
        clade_shift = paste(
            "A shifted response in a subtree or clade",
            "  Example: is there one branch of the tree with unusually strong KO responses?",
            sep = "\n"
        ),
        phylogenetic_mean = paste(
            "A phylogenetically corrected mean response",
            "  Example: what is the average response after accounting for shared ancestry?",
            sep = "\n"
        )
    )

    labels[.workflow_supported_goals(level)]
}

.workflow_variable_type_labels <- function() {
    c(
        categorical = paste(
            "Categorical variable, such as condition or treatment",
            "  Example: control versus warmed sediment",
            sep = "\n"
        ),
        continuous = paste(
            "Continuous variable, such as pH, salinity, or oxygen",
            "  Example: dissolved oxygen concentration across samples",
            sep = "\n"
        ),
        formula = paste(
            "Formula with additional covariates and one focal term",
            "  Example: test condition while adjusting for salinity and batch",
            sep = "\n"
        )
    )
}

.workflow_supported_goals <- function(level) {
    level <- .match_workflow_choice(level, .workflow_levels(), "level")

    switch(level,
        gene = "association",
        genome = c(
            "association",
            "group_coherence",
            "phylogenetic_signal",
            "clade_shift",
            "phylogenetic_mean"
        ),
        ko = c(
            "association",
            "heterogeneity",
            "group_difference",
            "group_coherence",
            "phylogenetic_signal",
            "clade_shift",
            "phylogenetic_mean"
        ),
        module = c(
            "coherent_effects",
            "total_activity",
            "heterogeneity",
            "group_difference",
            "phylogenetic_signal",
            "clade_shift",
            "phylogenetic_mean"
        ),
        pathway = c(
            "coherent_effects",
            "total_activity",
            "heterogeneity",
            "group_difference",
            "phylogenetic_signal",
            "clade_shift",
            "phylogenetic_mean"
        )
    )
}

.workflow_supports_phylo_count_model <- function(level, goal) {
    level %in% c("ko", "module", "pathway") &&
        goal %in% c("association", "heterogeneity", "group_difference", "total_activity")
}

.workflow_menu <- function(title, values, labels, menuFun = utils::menu) {
    choice <- menuFun(labels, title = title)

    if (!is.numeric(choice) || length(choice) != 1L || is.na(choice) || choice < 1L || choice > length(values)) {
        stop("Workflow selection cancelled.", call. = FALSE)
    }

    values[[choice]]
}

.workflow_prompt_logical <- function(title, trueLabel, falseLabel, menuFun = utils::menu) {
    .workflow_menu(
        title = title,
        values = list(FALSE, TRUE),
        labels = c(falseLabel, trueLabel),
        menuFun = menuFun
    )
}

.workflow_prompt_genome_offset <- function(menuFun = utils::menu) {
    .workflow_menu(
        title = paste(
            "Should the model include genome-abundance normalization",
            "through `genomeOffset`?"
        ),
        values = list(NULL, TRUE, FALSE),
        labels = c(
            paste(
                "Automatic: decide from the assays available in the object",
                "  Example: use this when you want the function to follow the experiment structure by default",
                sep = "\n"
            ),
            paste(
                "Yes: include genome-abundance normalization",
                "  Example: adjust KO or gene RNA for parent-genome abundance",
                sep = "\n"
            ),
            paste(
                "No: do not include genome-abundance normalization",
                "  Example: focus on RNA differences without abundance normalization",
                sep = "\n"
            )
        ),
        menuFun = menuFun
    )
}

.prompt_find_workflow <- function(
    level,
    goal,
    variableType,
    repeatedMeasures,
    genomeOffset,
    phylogeny,
    promptLevel,
    promptGoal,
    promptVariableType,
    promptRepeatedMeasures,
    promptGenomeOffset,
    promptPhylogeny,
    menuFun = utils::menu
) {
    if (!promptLevel && !is.null(level)) {
        level <- .match_workflow_choice(level, .workflow_levels(), "level")
    }

    if (isTRUE(promptLevel)) {
        level <- .workflow_menu(
            title = "What is your primary unit of interest?",
            values = as.list(.workflow_levels()),
            labels = unname(.workflow_level_labels()),
            menuFun = menuFun
        )
    }

    supported_goals <- .workflow_supported_goals(level)

    if (!promptGoal && !is.null(goal)) {
        goal <- .match_workflow_choice(goal, supported_goals, "goal")
    }

    if (isTRUE(promptGoal)) {
        goal <- .workflow_menu(
            title = "What biological question are you trying to answer?",
            values = as.list(supported_goals),
            labels = unname(.workflow_goal_labels(level)),
            menuFun = menuFun
        )
    }

    if (!promptVariableType && !is.null(variableType)) {
        variableType <- .match_workflow_choice(variableType, .workflow_variable_types(), "variableType")
    }

    if (isTRUE(promptVariableType)) {
        variableType <- .workflow_menu(
            title = "What type of predictor or model specification do you have?",
            values = as.list(.workflow_variable_types()),
            labels = unname(.workflow_variable_type_labels()),
            menuFun = menuFun
        )
    }

    if (isTRUE(promptRepeatedMeasures)) {
        repeatedMeasures <- .workflow_prompt_logical(
            title = "Do you need repeated-measures or blocked-sample support?",
            trueLabel = paste(
                "Yes: include a reminder about `sampleBlock`",
                "  Example: paired samples, time series, or several measurements from the same site",
                sep = "\n"
            ),
            falseLabel = paste(
                "No: no repeated-measures or sample block",
                "  Example: each sample is an independent observation",
                sep = "\n"
            ),
            menuFun = menuFun
        )
    }

    if (isTRUE(promptGenomeOffset)) {
        genomeOffset <- .workflow_prompt_genome_offset(menuFun = menuFun)
    }

    if (isTRUE(promptPhylogeny)) {
        phylogeny <- if (.workflow_supports_phylo_count_model(level, goal)) {
            .workflow_prompt_logical(
                title = "Should the count model account for phylogenetic non-independence among genomes?",
                trueLabel = paste(
                    "Yes: suggest `genomeCorrelation = \"brownian\"` when supported",
                    "  Example: closely related genomes are expected to have similar baseline activity or responses",
                    sep = "\n"
                ),
                falseLabel = paste(
                    "No: use independent genome effects",
                    "  Example: genome relatedness is not part of the current biological question",
                    sep = "\n"
                ),
                menuFun = menuFun
            )
        } else {
            FALSE
        }
    }

    list(
        level = level,
        goal = goal,
        variableType = variableType,
        repeatedMeasures = repeatedMeasures,
        genomeOffset = genomeOffset,
        phylogeny = phylogeny
    )
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

    if (!is.null(x$suggestedCode) && length(x$suggestedCode) > 0L) {
        cat("\nSuggested code:\n")
        for (one_line in x$suggestedCode) {
            cat(one_line, "\n", sep = "")
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
