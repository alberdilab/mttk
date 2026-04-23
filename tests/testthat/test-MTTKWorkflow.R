test_that("findWorkflow recommends gene-level association workflows", {
    rec <- findWorkflow(
        level = "gene",
        goal = "association",
        variableType = "categorical"
    )

    expect_s3_class(rec, "MTTKWorkflowRecommendation")
    expect_match(rec$question, "genes")
    expect_identical(as.character(rec$recommendations$workflow_function), "fitGeneModel")
    expect_identical(rec$suggestedArguments, "variable = \"...\"")
})

test_that("findWorkflow recommends genome group summaries", {
    rec <- findWorkflow(
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

test_that("findWorkflow recommends module coherent-effect workflows", {
    rec <- findWorkflow(
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

test_that("findWorkflow recommends KO phylogenetic workflows", {
    rec <- findWorkflow(
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

test_that("findWorkflow recommends pathway total-activity workflows", {
    rec <- findWorkflow(
        level = "pathway",
        goal = "total_activity",
        variableType = "categorical",
        genomeOffset = FALSE,
        phylogeny = TRUE
    )

    expect_identical(as.character(rec$recommendations$workflow_function), "fitPathwayMixedModel")
    expect_true(all(c("variable = \"...\"", "genomeOffset = FALSE", "genomeCorrelation = \"brownian\"") %in%
        rec$suggestedArguments))
})

test_that("findWorkflow rejects unsupported combinations", {
    expect_error(
        findWorkflow(level = "gene", goal = "heterogeneity"),
        "No mttk workflow is currently defined"
    )

    expect_error(
        findWorkflow(level = "module", goal = "group_coherence"),
        "No mttk workflow is currently defined"
    )
})
