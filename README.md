# MTTK

Metatranscriptomics toolkit for nested gene, genome, and community analysis.

## What MTTK is for

MTTK is designed for genome-resolved metatranscriptomics, where the data are
not only a gene-by-sample count matrix.

- genes are nested within genomes,
- genome abundance may be measured separately from gene RNA counts,
- genome relationships may be represented explicitly with a phylogenetic tree,
- genes may also belong to functional hierarchies such as KO, module, and
  pathway.

The package therefore focuses on:

- an explicit container that stores gene-level and genome-level measurements
  together,
- transparent link tables between genes, genomes, and annotations,
- aggregation and modeling workflows that can use those nested relationships.

## How MTTK Differs From edgeR, DESeq2, and limma-voom

`edgeR`, `DESeq2`, and `limma-voom` are established tools for differential
expression on feature-by-sample expression matrices. They are strong choices
when the main task is gene-wise or transcript-wise differential expression in a
mostly flat design.

MTTK is different in two main ways.

First, MTTK is built around a different data structure. It keeps:

- gene-level RNA data,
- genome-level data such as DNA abundance,
- explicit gene-to-genome nesting,
- optional functional hierarchies such as gene-to-KO, KO-to-module, and
  KO-to-pathway links, including many-to-many mappings

in one Bioconductor-friendly object.

Second, MTTK is aiming at questions that are awkward to express in a flat
RNA-seq workflow, for example:

- comparing activity while accounting for the parent genome of each gene,
- aggregating genes to KO or higher annotation levels before modeling,
- combining gene-level RNA with genome-level abundance offsets,
- moving between gene-, genome-, and function-level analyses in one framework.

The difference is also in the modeling target, not just in the container.

- `edgeR` typically fits one negative-binomial GLM per feature in a flat
  feature-by-sample design, with empirical-Bayes shrinkage for dispersion.
- `DESeq2` typically fits one negative-binomial GLM per feature in a flat
  feature-by-sample design, with shrinkage for dispersion and fold changes.
- `limma-voom` transforms counts to log-expression with precision weights, then
  fits weighted linear models per feature.
- MTTK's current function-level workflows start from a KO-level
  `glmmTMB` negative-binomial mixed model, then either:
  - summarize KO effect sizes upward to modules or pathways with a
    meta-analysis-like synthesis, or
  - re-aggregate counts and refit a mixed model at module or pathway level
    when the target is total functional activity.

So the current MTTK workflow answers a different question:

- standard RNA-seq workflows mostly ask whether each feature changes across
  samples,
- the current MTTK workflows can ask either whether KO-level effects assigned
  to a module or pathway are coherent, or whether total module/pathway activity
  shifts across samples.

Most modeling functions in MTTK support two user-facing interfaces:

- `variable = "condition"` for the simple one-variable case
- `formula = ~ condition + pH` plus `term = "condition"` when you want to
  adjust for additional sample covariates but still test one focal effect

## Choosing An Analysis

The higher-level question matters, because MTTK currently supports two
different module/pathway workflows.

- If the question is "which individual KOs are associated with the variable?",
  use `fitKOMixedModel()`.
- If the question is "does the same KO respond differently across genomes?",
  use `fitKORandomSlopeModel()` and then inspect the genome-specific KO effects
  with `koGenomeEffects()`.
- If the question is "does the KO response differ directly between two genome
  groups, such as Bacteria and Archaea?", use
  `fitKOGroupInteractionModel(..., group = "domain")`.
- If the question is "within a taxonomic group or clade, does this KO show a
  coherent genome-specific response?", fit KOs with random slopes first and
  then use `fitKOGenomeGroupMetaAnalysis(..., group = "domain")` or another
  grouping column from `genomeData(x)`.
- If the question is "does this KO respond differently between two clades or
  taxonomic groups?", first summarize KO genome effects with
  `fitKOGenomeGroupMetaAnalysis()` and then compare groups with
  `compareKOGenomeGroups()`.
- If the question is "do the KOs belonging to this module or pathway show a
  coherent association?", fit KOs first and then use
  `fitModuleMetaAnalysis()` or `fitPathwayMetaAnalysis()`.
  This works after either `fitKOMixedModel()` or `fitKORandomSlopeModel()`,
  because the meta-analysis layer summarizes the KO-level fixed-effect
  estimates from the fitted KO models.
  Then choose `membershipMode` according to the overlap problem:
  `duplicate` treats every annotated KO as fully belonging to every mapped set,
  `split` divides each KO contribution across its memberships to avoid
  over-counting broadly annotated KOs,
  and `exclusive` keeps only uniquely assigned KOs.
- If the question is "does the total RNA assigned to this module or pathway
  change?", use `fitModuleMixedModel()` or `fitPathwayMixedModel()`.
  Then choose `membershipMode = "duplicate"` if the target is total assigned
  activity, or `membershipMode = "exclusive"` if only uniquely assigned genes
  should contribute.
  `split` is intentionally not offered here because these workflows fit
  negative-binomial mixed models to counts, and splitting overlapping
  assignments would create fractional pseudo-counts rather than observed
  counts.
- If the question is "does the module or pathway response differ directly
  between two genome groups?", use `fitModuleGroupInteractionModel()` or
  `fitPathwayGroupInteractionModel()`.
- If the question is "does the same module or pathway respond differently
  across genomes?", use `fitModuleRandomSlopeModel()` or
  `fitPathwayRandomSlopeModel()` and then inspect the genome-specific
  conditional effects with `moduleGenomeEffects()` or
  `pathwayGenomeEffects()`.
- If the question is "which genomes respond?", use `fitGenomeModel()`.
- If the question is "do closely related genomes tend to respond similarly?",
  fit genomes first and then use `fitGenomePhylogeneticSignal()`.
- If the question is "is there a subtree of related genomes with a shifted
  response?", fit genomes first and then use `scanGenomeClades()`.
- If the question is "what is the average genome response after accounting for
  the phylogenetic covariance among genomes?", fit genomes first and then use
  `fitGenomePhylogeneticGLS()`.
- If the question is "do genomes in a clade or taxonomic group respond
  coherently?", fit genomes first and then use
  `fitGenomeGroupMetaAnalysis(..., group = "domain")` or another grouping
  column from `genomeData(x)`.
- If the question is "for a given KO, do genome-specific KO responses follow
  the phylogeny?", fit KOs with random slopes first and then use
  `fitKOPhylogeneticSignal()`.
- If the question is "is there a subtree where a KO responds in a distinct
  way?", fit KOs with random slopes first and then use `scanKOClades()`.
- If the question is "what is the average KO response after accounting for
  phylogenetic covariance among genomes?", fit KOs with random slopes first
  and then use `fitKOPhylogeneticGLS()`.
- If the question is "which individual genes change?", use `fitGeneModel()`.
- If the question is "do closely related genomes tend to show similar
  module- or pathway-level responses?", fit modules or pathways with random
  slopes first and then use `fitModulePhylogeneticSignal()` or
  `fitPathwayPhylogeneticSignal()`.
- If the question is "is there a responding subtree for a module or pathway?",
  fit modules or pathways with random slopes first and then use
  `scanModuleClades()` or `scanPathwayClades()`.
- If the question is "what is the average module or pathway response after
  accounting for phylogenetic covariance among genomes?", fit modules or
  pathways with random slopes first and then use
  `fitModulePhylogeneticGLS()` or `fitPathwayPhylogeneticGLS()`.

In practice, the KO meta-analysis route is the main higher-hierarchy option for
KEGG-style annotations, because KO membership in modules and pathways is often
many-to-many.

By default, the modeling workflows use a library-size offset and, when
`dna_genome_counts` is available in `genomeExperiment(x)`, a genome-abundance
offset. Use `genomeOffset = FALSE` when you want to keep library-size
normalization but disable genome-abundance normalization. In the KO/module/
pathway mixed-model workflows, the genome random effect is still included when
`genomeOffset = FALSE`.

The mode choice matters because it changes the exact question being answered:

- KO meta-analysis with `duplicate` asks whether the KOs annotated to a set,
  taken at face value, show a coherent effect.
- KO meta-analysis with `split` asks the same question while dividing each
  KO's contribution across its memberships, reducing the influence of KOs that
  participate in many sets.
- KO meta-analysis with `exclusive` asks whether uniquely assigned KOs alone
  support a set-level effect; this is the most specific and usually the most
  conservative option.
- Count reaggregation with `duplicate` asks whether total RNA assigned to a set
  changes, even if some genes are shared with other sets.
- Count reaggregation with `exclusive` asks whether the uniquely assigned part
  of a set changes, avoiding duplicated counts entirely but potentially
  discarding much of the data. There is no count-side `split` mode because the
  current reaggregation workflow uses negative-binomial count models, and
  splitting memberships would create fractional pseudo-counts.

So the current position of MTTK is:

- not a replacement for `edgeR`, `DESeq2`, or `limma-voom`,
- complementary infrastructure for nested metatranscriptomic data,
- a place to implement hierarchy-aware analyses that those packages were not
  primarily designed around.

## Current scope

The current version already supports:

- `MTTKExperiment` for explicit gene/genome data storage,
- aggregation to genomes and annotation-linked groups,
- KO-level mixed models with `glmmTMB`,
- KO-level random-slope mixed models with KO-by-genome effect extraction,
- module- and pathway-level KO meta-analysis built on KO effect sizes,
- module- and pathway-level mixed models for total functional activity,
- genome-level negative-binomial models,
- tree-based summaries of genome-level responses,
- phylogenetic GLS summaries of genome-level responses,
- tree-based summaries, clade scans, and phylogenetic GLS summaries for KO
  genome-specific responses,
- tree-based summaries, clade scans, and phylogenetic GLS summaries for
  module- and pathway-level genome-specific responses,
- genome-group meta-analysis across genome metadata columns such as domain or
  clade,
- gene-level negative-binomial models with optional parent-genome abundance
  offsets,
- optional KEGG module annotation retrieval through `KEGGREST` for joining
  module names and classes onto KO-level results.

Gene-wise differential expression with `edgeR`, `DESeq2`, or `limma-voom`
could still be used alongside MTTK, with MTTK handling the container,
alignment, and hierarchical summaries around those analyses.
