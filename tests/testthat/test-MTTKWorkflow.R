test_that("workflow builder recommends gene-level association workflows", {
    rec <- .build_workflow_recommendation(
        level = "gene",
        goal = "association",
        variableType = "categorical"
    )

    expect_s3_class(rec, "MTTKWorkflowRecommendation")
    expect_match(rec$question, "genes")
    expect_identical(as.character(rec$recommendations$workflow_function), "fitGeneModel")
    expect_identical(rec$suggestedArguments, "variable = \"...\"")
    expect_true(any(grepl("^gene_fit <- fitGeneModel\\($", rec$suggestedCode)))
    expect_true(any(grepl("variable = \"condition\"", rec$suggestedCode, fixed = TRUE)))
})

test_that("workflow builder recommends genome group summaries", {
    rec <- .build_workflow_recommendation(
        level = "genome",
        goal = "group_coherence",
        variableType = "continuous"
    )

    expect_identical(
        as.character(rec$recommendations$workflow_function),
        c("fitGenomeModel", "fitGenomeGroupMetaAnalysis")
    )
    expect_true(any(grepl("group =", rec$notes, fixed = TRUE)))
})

test_that("workflow builder recommends module coherent-effect workflows", {
    rec <- .build_workflow_recommendation(
        level = "module",
        goal = "coherent_effects",
        variableType = "continuous",
        phylogeny = TRUE
    )

    expect_identical(
        as.character(rec$recommendations$workflow_function),
        c("fitKOMixedModel", "fitModuleMetaAnalysis")
    )
    expect_true(any(grepl("membershipMode = \"split\"", rec$notes, fixed = TRUE)))
})

test_that("workflow builder recommends KO phylogenetic workflows", {
    rec <- .build_workflow_recommendation(
        level = "ko",
        goal = "phylogenetic_signal",
        variableType = "formula",
        repeatedMeasures = TRUE
    )

    expect_identical(
        as.character(rec$recommendations$workflow_function),
        c("fitKORandomSlopeModel", "fitKOPhylogeneticSignal")
    )
    expect_true(all(c("formula = ~ ... + ...", "term = \"...\"", "sampleBlock = \"...\"") %in%
        rec$suggestedArguments))
    expect_true(any(grepl("genomeTree", rec$notes, fixed = TRUE)))
})

test_that("workflow builder recommends pathway total-activity workflows", {
    rec <- .build_workflow_recommendation(
        level = "pathway",
        goal = "total_activity",
        variableType = "categorical",
        genomeOffset = FALSE,
        phylogeny = TRUE
    )

    expect_identical(as.character(rec$recommendations$workflow_function), "fitPathwayMixedModel")
    expect_true(all(c("variable = \"...\"", "genomeOffset = FALSE", "genomeCorrelation = \"brownian\"") %in%
        rec$suggestedArguments))
    expect_true(any(grepl("^pathway_fit <- fitPathwayMixedModel\\($", rec$suggestedCode)))
    expect_true(any(grepl("membershipMode = \"duplicate\"", rec$suggestedCode, fixed = TRUE)))
})

test_that("findWorkflow is interactive-only in non-interactive use", {
    expect_error(
        findWorkflow(),
        "interactive-only"
    )

    expect_error(
        findWorkflow(level = "ko"),
        "does not accept direct arguments"
    )
})

test_that("workflow recommendation print method shows suggested code", {
    rec <- .build_workflow_recommendation(
        level = "pathway",
        goal = "coherent_effects",
        variableType = "categorical"
    )

    output <- capture.output(print(rec))
    expect_true(any(grepl("Suggested code:", output, fixed = TRUE)))
    expect_true(any(grepl("fitPathwayMetaAnalysis", output, fixed = TRUE)))
})

test_that("interactive workflow helper walks through the decision tree", {
    answers <- c(3L, 2L, 3L, 2L, 1L, 2L)

    fake_menu <- function(choices, title) {
        out <- answers[[1]]
        answers <<- answers[-1]
        out
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
        promptPhylogeny = TRUE,
        menuFun = fake_menu
    )

    expect_identical(selected$level, "ko")
    expect_identical(selected$goal, "heterogeneity")
    expect_identical(selected$variableType, "formula")
    expect_true(selected$repeatedMeasures)
    expect_null(selected$genomeOffset)
    expect_true(selected$phylogeny)
})

test_that("workflow builder rejects unsupported combinations", {
    expect_error(
        .build_workflow_recommendation(level = "gene", goal = "heterogeneity", variableType = "categorical"),
        "No mttk workflow is currently defined"
    )

    expect_error(
        .build_workflow_recommendation(level = "module", goal = "group_coherence", variableType = "categorical"),
        "No mttk workflow is currently defined"
    )
})
