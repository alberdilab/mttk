#' Aggregate Gene RNA to KO-within-Genome Counts
#'
#' `aggregateToKOGenome()` collapses a gene-level assay to KO/genome pairs. The
#' returned rows represent KO-within-genome observations, which are the units
#' used by the first KO-level mixed-model workflow in MTTK.
#'
#' Row metadata always includes `ko_id`, `genome_id`, `n_genes_pair`, and
#' `n_links_pair`. Available `genomeData(x)` columns are appended when they can
#' be aligned unambiguously to the genome identifiers of the aggregated rows.
#'
#' @param x An `MTTKExperiment`.
#' @param assay A single gene-level assay name to aggregate.
#'
#' @return A `SummarizedExperiment` with one row per KO/genome pair and one
#'   column per sample.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' ko_genome <- aggregateToKOGenome(x)
#' ko_genome
#' SummarizedExperiment::rowData(ko_genome)[, c("ko_id", "genome_id")]
#'
#' @export
aggregateToKOGenome <- function(x, assay = "rna_gene_counts") {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    assay_name <- .normalize_analysis_assays(x, assay)
    if (length(assay_name) != 1L) {
        stop("'assay' must be a single gene-level assay name.", call. = FALSE)
    }

    .aggregate_rna_by_ko_genome(x, assay_name = assay_name)
}

#' Build a KO Mixed-Model Data Table
#'
#' `makeKOMixedModelData()` materializes the long-form observation table used by
#' [fitKOMixedModel()]. Each row corresponds to one KO/genome/sample
#' observation, with RNA counts, sample-level covariates, and optional
#' genome-abundance offsets aligned and ready for model fitting.
#'
#' The returned table is useful for inspecting the modeled data before fitting,
#' or for reusing the same KO/genome aggregation in custom workflows.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)` for the
#'   simple one-variable interface. The column must be numeric or a factor with
#'   at least two levels. For factors with three or more levels supply `term` to
#'   select which contrast to test. Supply exactly one of `variable` or `formula`.
#' @param formula Optional one-sided or two-sided fixed-effect formula for the
#'   sample-level covariates, for example `~ condition + pH` or
#'   `rna_count ~ condition + pH`. Offsets and genome-level random effects are
#'   added internally by MTTK. Supply exactly one of `variable` or `formula`.
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
#' @return An `S4Vectors::DataFrame` with one row per KO/genome/sample
#'   observation.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' model_data <- makeKOMixedModelData(x, variable = "condition")
#' model_data
#'
#' @export
makeKOMixedModelData <- function(
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

    model_spec <- .normalize_model_spec(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        referenceLevels = referenceLevels,
        allow_formula = TRUE,
        allow_random_slope = FALSE
    )

    out <- .make_feature_mixed_model_data(
        x = x,
        model_spec = model_spec,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = "duplicate",
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        path = "gene_to_ko",
        feature_label = "KO"
    )

    metadata_list <- S4Vectors::metadata(out)
    metadata_list$mttk_ko_model <- metadata_list$mttk_function_mixed_model
    S4Vectors::metadata(out) <- metadata_list
    out
}

#' Fit a KO-Level Mixed Model
#'
#' `fitKOMixedModel()` fits one negative-binomial mixed model per KO using
#' `glmmTMB`. RNA counts are first aggregated to KO-within-genome observations,
#' so the fitted model respects the nesting of genes within genomes.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`. The genome random
#' effect is part of this workflow regardless of the offset settings.
#' `genomeCorrelation = "independent"` uses the standard random intercept
#' `(1 | genome_id)`, while `genomeCorrelation = "brownian"` replaces it with a
#' Brownian phylogenetic covariance term derived from the genome tree.
#'
#' This workflow is intended for KO-level differential expression or association
#' analysis, where each KO can be observed in multiple genomes and the genome
#' term should be modeled explicitly as a random intercept.
#'
#' @param x An `MTTKExperiment`.
#' @param variable A single sample-level column name from `colData(x)` for the
#'   simple one-variable interface. The column must be numeric or a factor with
#'   at least two levels. For factors with three or more levels supply `term` to
#'   select which contrast to test. Supply exactly one of `variable` or `formula`.
#' @param formula Optional one-sided or two-sided fixed-effect formula for the
#'   sample-level covariates, for example `~ condition + salinity` or
#'   `rna_count ~ condition + salinity`. Offsets and genome-level random
#'   effects are added internally by MTTK. Supply exactly one of `variable` or
#'   `formula`.
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
#' @param genomeCorrelation Genome-level correlation structure for the random
#'   intercept. Use `"independent"` for the standard genome random intercept, or
#'   `"brownian"` to model genome-level covariance from a phylogenetic tree via
#'   `glmmTMB::propto()`.
#' @param phylogeny Optional `ape::phylo` object used when
#'   `genomeCorrelation = "brownian"`. If `NULL`, MTTK uses `genomeTree(x)`.
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
#' @param BPPARAM A `BiocParallelParam` instance controlling how models are
#'   fitted across KOs. Defaults to `BiocParallel::SerialParam()` for
#'   single-core execution. Pass e.g. `BiocParallel::MulticoreParam(4)` to fit
#'   KO models in parallel.
#'
#' @return An `MTTKFit` with one row per KO.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitKOMixedModel(x, variable = "condition")
#'     fit
#'     significantResults(fit)
#' }
#'
#' @export
fitKOMixedModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
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
    .fit_feature_mixed_model(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = "duplicate",
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        keepFits = keepFits,
        path = "gene_to_ko",
        feature_label = "KO",
        fit_model_name = "ko_mixed_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a KO-Level Genome-Group Interaction Model
#'
#' `fitKOGroupInteractionModel()` fits one negative-binomial mixed model per KO
#' using `glmmTMB`, with a genome random intercept and a fixed interaction
#' between the selected sample-level effect and a two-level genome grouping
#' variable from `genomeData(x)`.
#'
#' RNA counts are first aggregated to KO-within-genome observations. The fitted
#' interaction term answers the direct question: does the KO response differ
#' between two genome groups, such as domains or selected clades, after
#' adjusting for the chosen sample-level covariates and offsets?
#'
#' @inheritParams fitKOMixedModel
#' @param group A single column name from `genomeData(x)` defining the genome
#'   groups to compare.
#' @param groupLevels Optional character vector of length 2 giving the genome
#'   group reference and contrast levels. If `NULL`, the grouping column must
#'   contain exactly two usable groups.
#'
#' @return An `MTTKFit` with one row per KO.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     fit <- fitKOGroupInteractionModel(
#'         x,
#'         variable = "condition",
#'         group = "domain"
#'     )
#'     fit
#' }
#'
#' @export
fitKOGroupInteractionModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    group,
    groupLevels = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
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
        membershipMode = "duplicate",
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        keepFits = keepFits,
        path = "gene_to_ko",
        feature_label = "KO",
        fit_model_name = "ko_group_interaction_model",
        BPPARAM = BPPARAM
    )
}

#' Fit a KO-Level Random-Slope Mixed Model
#'
#' `fitKORandomSlopeModel()` fits one negative-binomial mixed model per KO using
#' `glmmTMB`, with both a random intercept and a random slope for the selected
#' sample-level variable across genomes.
#'
#' RNA counts are first aggregated to KO-within-genome observations. The fitted
#' model then estimates an overall KO-level effect together with genome-specific
#' deviations around that effect. This workflow is intended for questions such
#' as whether the same KO responds differently across genomes, and whether those
#' genome-specific KO responses appear similar or heterogeneous.
#'
#' `libraryOffset` controls whether the model includes
#' `offset(log(lib_size))`. `genomeOffset` controls whether the model includes
#' `offset(log(genome_abundance + offsetPseudocount))`.
#' `genomeCorrelation = "independent"` uses the standard random intercept and
#' random slope across genomes, while `genomeCorrelation = "brownian"` replaces
#' both with a joint Brownian phylogenetic intercept-slope covariance structure
#' derived from the genome tree.
#'
#' Genome-specific KO coefficients can be extracted with [koGenomeEffects()].
#'
#' @inheritParams fitKOMixedModel
#' @param randomSlope Optional sample-level variable whose genome-specific slope
#'   should be modeled when `formula` is used. If omitted, MTTK uses the tested
#'   variable when it can be inferred unambiguously.
#'
#' @return An `MTTKFit` with one row per KO. The returned fit also stores
#'   KO-by-genome conditional effects that can be extracted with
#'   [koGenomeEffects()].
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     fit
#'     koGenomeEffects(fit)
#' }
#'
#' @export
fitKORandomSlopeModel <- function(
    x,
    variable = NULL,
    formula = NULL,
    term = NULL,
    assay = "rna_gene_counts",
    libraryOffset = TRUE,
    genomeOffset = NULL,
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
    .fit_feature_random_slope_model(
        x = x,
        variable = variable,
        formula = formula,
        term = term,
        assay = assay,
        libraryOffset = libraryOffset,
        genomeOffset = genomeOffset,
        membershipMode = "duplicate",
        libSize = libSize,
        sampleBlock = sampleBlock,
        genomeCorrelation = genomeCorrelation,
        phylogeny = phylogeny,
        genomeAssay = genomeAssay,
        offsetPseudocount = offsetPseudocount,
        referenceLevels = referenceLevels,
        randomSlope = randomSlope,
        keepFits = keepFits,
        path = "gene_to_ko",
        feature_label = "KO",
        fit_model_name = "ko_random_slope_model",
        BPPARAM = BPPARAM
    )
}
