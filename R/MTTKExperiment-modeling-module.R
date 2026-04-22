#' Aggregate Gene RNA to Module-within-Genome Counts
#'
#' `aggregateToModuleGenome()` collapses a gene-level assay to module/genome
#' pairs by following the functional path `gene_to_ko -> ko_to_module`. The
#' returned rows represent module-within-genome observations that can be used in
#' module-level mixed models.
#'
#' Because KO-to-module membership can be many-to-many, the same gene can
#' contribute to multiple modules. `membershipMode = "duplicate"` carries the
#' full count into each mapped module/genome observation and is appropriate for
#' the question "how much RNA is assigned to this module?". In contrast,
#' `membershipMode = "exclusive"` keeps only genes whose final module mapping is
#' unique, which is more conservative and avoids duplicated counts entirely.
#'
#' @param x An `MTTKExperiment`.
#' @param assay A single gene-level assay name to aggregate.
#' @param membershipMode How to handle many-to-many set memberships.
#'   `"duplicate"` carries a gene's RNA into every mapped module, while
#'   `"exclusive"` keeps only uniquely mapped genes.
#'
#' @return A `SummarizedExperiment` with one row per module/genome pair and one
#'   column per sample.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' module_genome <- aggregateToModuleGenome(x)
#' module_genome
#' SummarizedExperiment::rowData(module_genome)[, c("module_id", "genome_id")]
#'
#' @export
aggregateToModuleGenome <- function(
    x,
    assay = "rna_gene_counts",
    membershipMode = c("duplicate", "exclusive")
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    assay_name <- .normalize_analysis_assays(x, assay)
    if (length(assay_name) != 1L) {
        stop("'assay' must be a single gene-level assay name.", call. = FALSE)
    }

    .aggregate_rna_by_feature_genome(
        x = x,
        assay_name = assay_name,
        path = c("gene_to_ko", "ko_to_module"),
        feature_label = "module",
        membership_mode = match.arg(membershipMode)
    )
}

#' Build a Module Mixed-Model Data Table
#'
#' `makeModuleMixedModelData()` materializes the long-form observation table
#' used by [fitModuleMixedModel()]. Each row corresponds to one
#' module/genome/sample observation, with RNA counts, sample-level covariates,
#' and optional genome-abundance offsets aligned and ready for model fitting.
#'
#' @inheritParams makeKOMixedModelData
#' @inheritParams aggregateToModuleGenome
#'
#' @return An `S4Vectors::DataFrame` with one row per module/genome/sample
#'   observation.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' model_data <- makeModuleMixedModelData(x, variable = "condition")
#' model_data
#'
#' @export
makeModuleMixedModelData <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
    membershipMode = c("duplicate", "exclusive"),
    libSize = NULL,
    sampleBlock = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL
) {
    membershipMode <- match.arg(membershipMode)
    model_spec <- .normalize_model_spec(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        referenceLevels = referenceLevels,
        allow_formula = TRUE,
        allow_random_slope = FALSE
    )
    .make_feature_mixed_model_data(
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
        path = c("gene_to_ko", "ko_to_module"),
        feature_label = "module"
    )
}

#' Fit a Module-Level Mixed Model
#'
#' `fitModuleMixedModel()` fits one negative-binomial mixed model per module
#' using `glmmTMB`. Gene-level RNA counts are first aggregated to
#' module-within-genome observations by following `gene_to_ko -> ko_to_module`,
#' then each module is modeled across genomes and samples.
#'
#' This workflow is intended for the total-activity question: does the summed
#' RNA assigned to a module shift across samples after accounting for genome
#' structure and optional genome abundance? For the different question of
#' whether the KO-level effects assigned to a module show a coherent direction
#' of change, use [fitModuleMetaAnalysis()] instead.
#'
#' `membershipMode = "duplicate"` answers that question in terms of all RNA
#' assigned to the module, even if the same KO contributes to multiple modules.
#' `membershipMode = "exclusive"` instead restricts the analysis to uniquely
#' assigned memberships, which avoids duplicated counts but may remove many
#' genes when module overlap is dense.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`. The genome random
#' effect `(1 | genome_id)` is part of this workflow regardless of the offset
#' settings.
#'
#' @inheritParams fitKOMixedModel
#' @inheritParams aggregateToModuleGenome
#'
#' @return An `MTTKFit` with one row per module.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitModuleMixedModel(x, variable = "condition")
#'     fit
#' }
#'
#' @export
fitModuleMixedModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
    membershipMode = c("duplicate", "exclusive"),
    libSize = NULL,
    sampleBlock = NULL,
    genomeCorrelation = c("independent", "brownian"),
    phylogeny = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL,
    keepFits = FALSE,
    BPPARAM = BiocParallel::SerialParam()
) {
    membershipMode <- match.arg(membershipMode)
    .fit_feature_mixed_model(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        keepFits = keepFits,
        path = c("gene_to_ko", "ko_to_module"),
        feature_label = "module",
        fit_model_name = "module_mixed_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a Module-Level Genome-Group Interaction Model
#'
#' `fitModuleGroupInteractionModel()` fits one negative-binomial mixed model per
#' module using `glmmTMB`, with a genome random intercept and a fixed
#' interaction between the selected sample-level effect and a two-level genome
#' grouping variable from `genomeData(x)`.
#'
#' This workflow is intended for the direct question: does the module response
#' differ between two genome groups after accounting for genome nesting,
#' optional genome abundance, and any covariates included in the fixed-effects
#' formula?
#'
#' @inheritParams fitModuleMixedModel
#' @param group A single column name from `genomeData(x)` defining the genome
#'   groups to compare.
#' @param groupLevels Optional character vector of length 2 giving the genome
#'   group reference and contrast levels. If `NULL`, the grouping column must
#'   contain exactly two usable groups.
#'
#' @return An `MTTKFit` with one row per module.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     fit <- fitModuleGroupInteractionModel(
#'         x,
#'         variable = "condition",
#'         group = "domain"
#'     )
#'     fit
#' }
#'
#' @export
fitModuleGroupInteractionModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    group,
    groupLevels = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
    membershipMode = c("duplicate", "exclusive"),
    libSize = NULL,
    sampleBlock = NULL,
    genomeCorrelation = c("independent", "brownian"),
    phylogeny = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL,
    keepFits = FALSE,
    BPPARAM = BiocParallel::SerialParam()
) {
    membershipMode <- match.arg(membershipMode)
    .fit_feature_group_interaction_model(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        group = group,
        groupLevels = groupLevels,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        keepFits = keepFits,
        path = c("gene_to_ko", "ko_to_module"),
        feature_label = "module",
        fit_model_name = "module_group_interaction_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a Module-Level Random-Slope Mixed Model
#'
#' `fitModuleRandomSlopeModel()` fits one negative-binomial mixed model per
#' module using `glmmTMB`, with both a genome random intercept and a
#' genome-specific random slope for the selected sample variable.
#'
#' Gene-level RNA counts are first aggregated to module-within-genome
#' observations by following `gene_to_ko -> ko_to_module`. The fitted model then
#' estimates an overall module-level effect together with genome-specific
#' deviations around that module effect.
#'
#' This workflow is intended for the question: does the same module show similar
#' or heterogeneous condition-associated responses across genomes? For the
#' higher-level question of whether KO-level module members show a coherent
#' direction of change without re-aggregating counts, use
#' [fitModuleMetaAnalysis()] instead.
#'
#' `membershipMode = "duplicate"` answers that question in terms of all RNA
#' assigned to the module, even if the same KO contributes to multiple modules.
#' `membershipMode = "exclusive"` instead restricts the analysis to uniquely
#' assigned memberships.
#'
#' `genomeCorrelation = "independent"` uses the standard random intercept and
#' random slope across genomes, while `genomeCorrelation = "brownian"` replaces
#' both with a joint Brownian phylogenetic intercept-slope covariance structure
#' derived from the genome tree.
#'
#' Genome-specific conditional module effects can be extracted with
#' [moduleGenomeEffects()].
#'
#' @inheritParams fitModuleMixedModel
#' @param randomSlope Optional sample-level variable whose genome-specific slope
#'   should be estimated. When `variable` is supplied, MTTK uses that variable.
#'   When `formula` is supplied, `randomSlope` can be omitted only if MTTK can
#'   infer a single tested variable unambiguously.
#'
#' @return An `MTTKFit` with one row per module. The returned fit also stores
#'   module-by-genome conditional effects that can be extracted with
#'   [moduleGenomeEffects()].
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitModuleRandomSlopeModel(x, variable = "condition")
#'     fit
#'     moduleGenomeEffects(fit)
#' }
#'
#' @export
fitModuleRandomSlopeModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
    membershipMode = c("duplicate", "exclusive"),
    libSize = NULL,
    sampleBlock = NULL,
    genomeCorrelation = c("independent", "brownian"),
    phylogeny = NULL,
    genomeAssay = "dna_genome_counts",
    offsetPseudocount = 1,
    referenceLevels = NULL,
    randomSlope = NULL,
    keepFits = FALSE,
    BPPARAM = BiocParallel::SerialParam()
) {
    membershipMode <- match.arg(membershipMode)
    .fit_feature_random_slope_model(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = membershipMode,
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        randomSlope = randomSlope,
        keepFits = keepFits,
        path = c("gene_to_ko", "ko_to_module"),
        feature_label = "module",
        fit_model_name = "module_random_slope_model",
        BPPARAM = BPPARAM
    )
}
