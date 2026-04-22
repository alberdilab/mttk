#' Aggregate Gene RNA to Pathway-within-Genome Counts
#'
#' `aggregateToPathwayGenome()` collapses a gene-level assay to pathway/genome
#' pairs by following the functional path `gene_to_ko -> ko_to_pathway`. The
#' returned rows represent pathway-within-genome observations that can be used
#' in pathway-level mixed models.
#'
#' Because KO-to-pathway membership can be many-to-many, the same gene can
#' contribute to multiple pathways. `membershipMode = "duplicate"` carries the
#' full count into each mapped pathway/genome observation and is appropriate
#' for the question "how much RNA is assigned to this pathway?". In contrast,
#' `membershipMode = "exclusive"` keeps only genes whose final pathway mapping
#' is unique, which is more conservative and avoids duplicated counts entirely.
#'
#' @param x An `MTTKExperiment`.
#' @param assay A single gene-level assay name to aggregate.
#' @param membershipMode How to handle many-to-many set memberships.
#'   `"duplicate"` carries a gene's RNA into every mapped pathway, while
#'   `"exclusive"` keeps only uniquely mapped genes.
#'
#' @return A `SummarizedExperiment` with one row per pathway/genome pair and one
#'   column per sample.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' pathway_genome <- aggregateToPathwayGenome(x)
#' pathway_genome
#' SummarizedExperiment::rowData(pathway_genome)[, c("pathway_id", "genome_id")]
#'
#' @export
aggregateToPathwayGenome <- function(
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
        path = c("gene_to_ko", "ko_to_pathway"),
        feature_label = "pathway",
        membership_mode = match.arg(membershipMode)
    )
}

#' Build a Pathway Mixed-Model Data Table
#'
#' `makePathwayMixedModelData()` materializes the long-form observation table
#' used by [fitPathwayMixedModel()]. Each row corresponds to one
#' pathway/genome/sample observation, with RNA counts, sample-level covariates,
#' and optional genome-abundance offsets aligned and ready for model fitting.
#'
#' @inheritParams makeKOMixedModelData
#' @inheritParams aggregateToPathwayGenome
#'
#' @return An `S4Vectors::DataFrame` with one row per pathway/genome/sample
#'   observation.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' model_data <- makePathwayMixedModelData(x, variable = "condition")
#' model_data
#'
#' @export
makePathwayMixedModelData <- function(
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
        path = c("gene_to_ko", "ko_to_pathway"),
        feature_label = "pathway"
    )
}

#' Fit a Pathway-Level Mixed Model
#'
#' `fitPathwayMixedModel()` fits one negative-binomial mixed model per pathway
#' using `glmmTMB`. Gene-level RNA counts are first aggregated to
#' pathway-within-genome observations by following `gene_to_ko -> ko_to_pathway`,
#' then each pathway is modeled across genomes and samples.
#'
#' This workflow is intended for the total-activity question: does the summed
#' RNA assigned to a pathway shift across samples after accounting for genome
#' structure and optional genome abundance? For the different question of
#' whether the KO-level effects assigned to a pathway show a coherent direction
#' of change, use [fitPathwayMetaAnalysis()] instead.
#'
#' `membershipMode = "duplicate"` answers that question in terms of all RNA
#' assigned to the pathway, even if the same KO contributes to multiple
#' pathways. `membershipMode = "exclusive"` instead restricts the analysis to
#' uniquely assigned memberships, which avoids duplicated counts but may remove
#' many genes when pathway overlap is dense.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`. The genome random
#' effect `(1 | genome_id)` is part of this workflow regardless of the offset
#' settings.
#'
#' @inheritParams fitKOMixedModel
#' @inheritParams aggregateToPathwayGenome
#'
#' @return An `MTTKFit` with one row per pathway.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitPathwayMixedModel(x, variable = "condition")
#'     fit
#' }
#'
#' @export
fitPathwayMixedModel <- function(
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
        path = c("gene_to_ko", "ko_to_pathway"),
        feature_label = "pathway",
        fit_model_name = "pathway_mixed_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a Pathway-Level Genome-Group Interaction Model
#'
#' `fitPathwayGroupInteractionModel()` fits one negative-binomial mixed model
#' per pathway using `glmmTMB`, with a genome random intercept and a fixed
#' interaction between the selected sample-level effect and a two-level genome
#' grouping variable from `genomeData(x)`.
#'
#' This workflow is intended for the direct question: does the pathway response
#' differ between two genome groups after accounting for genome nesting,
#' optional genome abundance, and any covariates included in the fixed-effects
#' formula?
#'
#' @inheritParams fitPathwayMixedModel
#' @param group A single column name from `genomeData(x)` defining the genome
#'   groups to compare.
#' @param groupLevels Optional character vector of length 2 giving the genome
#'   group reference and contrast levels. If `NULL`, the grouping column must
#'   contain exactly two usable groups.
#'
#' @return An `MTTKFit` with one row per pathway.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     fit <- fitPathwayGroupInteractionModel(
#'         x,
#'         variable = "condition",
#'         group = "domain"
#'     )
#'     fit
#' }
#'
#' @export
fitPathwayGroupInteractionModel <- function(
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
        path = c("gene_to_ko", "ko_to_pathway"),
        feature_label = "pathway",
        fit_model_name = "pathway_group_interaction_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a Pathway-Level Random-Slope Mixed Model
#'
#' `fitPathwayRandomSlopeModel()` fits one negative-binomial mixed model per
#' pathway using `glmmTMB`, with both a genome random intercept and a
#' genome-specific random slope for the selected sample variable.
#'
#' Gene-level RNA counts are first aggregated to pathway-within-genome
#' observations by following `gene_to_ko -> ko_to_pathway`. The fitted model
#' then estimates an overall pathway-level effect together with genome-specific
#' deviations around that pathway effect.
#'
#' This workflow is intended for the question: does the same pathway show
#' similar or heterogeneous condition-associated responses across genomes? For
#' the different higher-level question of whether KO-level pathway members show
#' a coherent direction of change without re-aggregating counts, use
#' [fitPathwayMetaAnalysis()] instead.
#'
#' `membershipMode = "duplicate"` answers that question in terms of all RNA
#' assigned to the pathway, even if the same KO contributes to multiple
#' pathways. `membershipMode = "exclusive"` instead restricts the analysis to
#' uniquely assigned memberships.
#'
#' `genomeCorrelation = "independent"` uses the standard random intercept and
#' random slope across genomes, while `genomeCorrelation = "brownian"` replaces
#' both with a joint Brownian phylogenetic intercept-slope covariance structure
#' derived from the genome tree.
#'
#' Genome-specific conditional pathway effects can be extracted with
#' [pathwayGenomeEffects()].
#'
#' @inheritParams fitPathwayMixedModel
#' @param randomSlope Optional sample-level variable whose genome-specific slope
#'   should be estimated. When `variable` is supplied, MTTK uses that variable.
#'   When `formula` is supplied, `randomSlope` can be omitted only if MTTK can
#'   infer a single tested variable unambiguously.
#'
#' @return An `MTTKFit` with one row per pathway. The returned fit also stores
#'   pathway-by-genome conditional effects that can be extracted with
#'   [pathwayGenomeEffects()].
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitPathwayRandomSlopeModel(x, variable = "condition")
#'     fit
#'     pathwayGenomeEffects(fit)
#' }
#'
#' @export
fitPathwayRandomSlopeModel <- function(
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
        path = c("gene_to_ko", "ko_to_pathway"),
        feature_label = "pathway",
        fit_model_name = "pathway_random_slope_model",
        BPPARAM = BPPARAM
    )
}
