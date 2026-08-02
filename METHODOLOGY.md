# Methodology

## Comparative analysis of mammalian gut microbiome functions associated with species-level cancer mortality

This document expands the project's 9-phase workflow (Phase 0–8) into a detailed, step-by-step methodology, with current best practices and their sources cited at each step. For each phase, it also covers the main alternative approaches that were considered and why they were set aside in favor of the chosen approach, given this project's specific constraints (small species-level *n*, multi-order comparative scope, and existing Python/R tooling). It is meant to function as the methods section of the project README: someone unfamiliar with the project should be able to follow it and understand not just *what* to do at each step, but *why*, what else was on the table, and what tool/paper each choice is based on.

Bracketed numbers `[n]` refer to the numbered reference list at the end of this document.

---

## Phase 0 — Scope the question

**Objective:** Fix one specific, falsifiable hypothesis before any data is touched.

### Steps

1. **Write the hypothesis as a single testable sentence**, not a topic. "Gut microbiome relates to cancer" is a topic; "species with a higher relative abundance of genes for secondary bile acid 7α-dehydroxylation have higher cancer mortality risk, after correcting for phylogeny and body mass" is a hypothesis. The distinction matters because each candidate endpoint requires a different data product downstream:
   - **Taxonomic composition** endpoint → needs only relative-abundance tables (16S or shotgun taxonomic profiles).
   - **Functional/pathway content** endpoint → needs gene-level annotation of assembled genomes (KEGG/CAZy), one level of processing beyond taxonomy.
   - **Predicted metabolite output** endpoint → needs genome-scale metabolic models and flux balance analysis (FBA), the most processing-intensive and most mechanistically informative option.

   Given the project's interest in specific oncogenic/protective metabolites (bile acids, H₂S, colibactin, SCFAs), the functional/metabolite-level endpoint is the appropriate target, but it should be stated explicitly rather than left implicit, because it determines every downstream tool choice.

2. **Decide the covariates that must be in every model up front, not post hoc.** The two cancer-mortality source papers already establish that litter size [2] and body mass/lifespan (tested and *not* found significant, i.e. Peto's paradox holds) are relevant life-history covariates. Deciding this before seeing the microbiome data prevents the covariate list from being quietly adjusted after seeing which model "worked" — a form of researcher degrees of freedom that inflates false-positive rates.

3. **Write this down before pulling data.** Formal preregistration — publicly time-stamping the hypothesis, covariates, and planned statistical test before the analysis is run — is standard practice for reducing exactly this kind of degrees-of-freedom problem in hypothesis-driven biology [6]. For a project this size, a lightweight version (a dated section in this repository, committed before Phase 1 begins) captures most of the benefit without the overhead of a formal registry.

### Best-practice rationale
Preregistration works by fixing the analysis plan *before* the outcome is known, which removes the ability (even unintentional) to choose the hypothesis, covariates, or test that produces the most favorable-looking result after the fact [6]. This is the single cheapest thing to do in the whole pipeline and the easiest to skip, which is why it is Phase 0.

### Alternatives considered
- **Purely exploratory / data-driven scan, with no fixed hypothesis.** Instead of committing to a specific endpoint and covariate set, one could compute correlations between every microbiome feature and cancer mortality and see what turns up. This was not chosen because with only ~55 species and potentially hundreds of candidate features, an unconstrained scan is exactly the maximal researcher-degrees-of-freedom scenario — with that many comparisons and that few data points, some features will correlate with cancer mortality by chance alone, and without a pre-fixed hypothesis there's no principled way to tell a real signal from noise dressed up as a discovery.
- **Full formal preregistration via a public registry (e.g., OSF).** This is the more rigorous version of Phase 0's lightweight approach and is the better choice once the project is closer to a publishable analysis. It wasn't chosen for the current build-out stage because it adds process overhead (registry submission, harder-to-revise commitments) that isn't justified while the pipeline itself is still being assembled and debugged; a dated, version-controlled statement in this repository captures most of the same benefit and can be upgraded to a formal OSF preregistration before the Phase 5 analysis is run for real.

### Deliverable
A one-paragraph, dated hypothesis statement + covariate list, committed to the repository before Phase 1 work starts.

---

## Phase 1 — Data acquisition and species overlap check

**Objective:** Determine the real, usable sample size before any bioinformatics work, since this is the single biggest feasibility gate on the project.

### Steps

1. **Pull each source dataset's species list independently.** For this project: cancer mortality from Vincze et al. 2022 [1] (191 species, via the authors' public GitHub repository) and Boddy et al. 2020 [2] (37 species, via journal supplementary Table S2); microbiome/metabolome data from Milani et al. 2020 [3] (77 species), Youngblut et al. 2020 [4] (180 species across 5 vertebrate classes, 103 of which are mammals), and Gregor et al. 2022 [5] (25 species).

2. **Standardize taxonomic names before joining anything.** This is the step most likely to silently corrupt a merge. Species names arrive in inconsistent formats across papers — underscore- vs. space-separated (`Equus_africanus` vs. `Equus africanus`), trinomial subspecies vs. binomial species (`Equus africanus asinus` vs. `Equus africanus`), and outright taxonomic synonyms (the same organism under two different scientific names due to a taxonomic revision). A purely string-based normalizer (lowercase, strip underscores, collapse trinomials to binomials) catches the first two problems but not synonyms. Catching synonyms requires cross-referencing against a taxonomic backbone — GBIF's Backbone Taxonomy or NCBI Taxonomy are the two most widely used [7]. The R package `taxize` wraps both of these (and several other) taxonomic web services in a single interface specifically for this kind of name-resolution task in ecological/biological datasets [7].

3. **Compute the overlap as an actual intersection, species by species, and log it.** Build one row per species with a column per source dataset, then define "usable" as species present in at least one cancer-mortality source *and* at least one microbiome/metabolome source. Do not estimate this from the reported *n*s in each paper's abstract — actually joining on species names is the only way to find out, since the overlap is frequently much smaller than either total looks. (For this project, real joins found 55 overlapping species out of a combined pool of ~400 unique names — about 14–28% of any individual source's species list, depending on which source you start from.)

4. **Treat the resulting number as a go/no-go gate.** If the overlap is too thin to support the planned statistical model (see Phase 8 for what "too thin" means quantitatively), the options are: (a) generate new sequencing data for specific missing species, (b) narrow scope to the best-covered mammalian orders (e.g., Carnivora and Primates are heavily represented across all five sources here), or (c) retrieve additional public datasets to widen coverage. This decision should happen before Phase 2, not after a failed Phase 5 model.

### Best-practice rationale
Cross-dataset species matching is a known, well-documented failure point in comparative biology specifically because taxonomic nomenclature is not perfectly standardized across independent data sources, and mismatches fail silently (a misspelled or synonym name simply drops out of the join instead of throwing an error) [7]. Explicitly normalizing names and logging the before/after counts at each step (as this project's `phase1_data_ingestion.py` does) is the only way to catch this before it propagates into a biased or undersized downstream analysis.

### Alternatives considered
- **Resolve every species name through a full taxonomic-backbone service (GBIF/NCBI via `taxize`) as the first pass, rather than a lightweight string normalizer.** Full backbone resolution is more thorough — it catches synonyms immediately — but it requires a network call per species, adds an external dependency, and is slower to iterate on during pipeline development. The lightweight string normalizer is used first because it's fast, offline, and catches the majority of real-world mismatches (formatting and subspecies-vs-species differences); `taxize`-based resolution is reserved for the smaller number of names that still fail to match after normalization, and Phase 2's tree-pruning step doubles as a second, independent check for whatever slips through both passes.
- **Include additional cancer-mortality sources beyond Vincze et al. and Boddy et al.** — for example, zoo/aquarium necropsy registries such as Species360's ZIMS database. These likely contain more species and richer necropsy detail, but access requires institutional data-use agreements and the records aren't published with the standardized species-level summary statistics (CMR, neoplasia prevalence) the two chosen papers already provide freely. Vincze et al. and Boddy et al. were chosen as the two largest publicly accessible, already-published, species-level cancer datasets that can be joined without a data-access negotiation.
- **Use only the single largest microbiome dataset (Youngblut et al., 103 mammal species) rather than combining three microbiome sources.** This would simplify the merge logic considerably. It wasn't chosen because restricting to one source caps the usable overlap at whatever intersects Youngblut's species list alone, which is smaller than the 55-species overlap achieved by combining all three microbiome sources — and sample size is the primary constraint on this project's statistical power (see Phase 8), so the extra merge complexity is worth the larger *n*.

### Deliverable
A species-by-source matrix and an overlap-only species list (already produced for this project — see `output/species_overlap_report.csv` and `output/overlap_species_only.csv`), plus a short written go/no-go decision based on the resulting sample size.

---

## Phase 2 — Build a species-level phylogenetic tree

**Objective:** Obtain a phylogeny containing exactly the overlap species from Phase 1, since every downstream statistical step (Phase 5) requires it.

### Steps

1. **Pull a tree for the full overlap species list from a maintained divergence-time database.** TimeTree.org aggregates published molecular divergence-time estimates into a single queryable resource and is the standard public source for this [8]. For mammal-specific work, the VertLife project's mammal supertree (Upham et al. 2019) is the more commonly used source in recent mammalian comparative-cancer literature — it provides not a single tree but a *credible set* of trees capturing topological and divergence-time uncertainty across ~6,000 living mammal species, built from a combination of DNA-sequence data (for the ~5,000 species with usable sequence data) and taxonomy-informed placement (for the remainder) [9].

2. **Do not use a taxonomy-only "birth-death polytomy resolver" tree for trait-based analyses if you can avoid it.** Trees that place unsampled species onto a backbone phylogeny using taxonomic rank alone (rather than sequence data) position those species at effectively random locations with respect to any biological trait — which breaks the very phylogenetic signal PGLS is trying to model, and biases the results in ways that are hard to detect after the fact [9]. Upham et al.'s credible-tree-set approach is preferred specifically because it is explicit about which species are sequence-placed versus taxonomy-placed, so this can be checked.

3. **Prune the tree to the exact overlap species list, both directions.** Any tree tip not in the overlap-species dataset should be dropped (`drop.tip()` in the R `ape` package is the standard tool for this), and — just as important — any overlap species *not* found in the tree should be flagged and resolved (usually a naming mismatch caught late) rather than silently dropped from the analysis.

4. **Use this same pruning step as a second, independent taxonomic QC pass.** Species names that fail to match the tree's tip labels despite matching across the Phase 1 datasets are a strong signal of a remaining synonym problem Phase 1's string-based normalizer missed — this is a good, cheap place to catch what the simpler QC step in Phase 1 couldn't.

### Best-practice rationale
PGLS requires a fully bifurcating (or at least branch-length-resolved) tree over exactly the species being modeled; unresolved nodes (polytomies) or a tree topology inferred independently of sequence data for some species will bias the phylogenetic correction PGLS depends on [9]. Building this tree early — before committing to bioinformatics work in Phase 3 — means naming inconsistencies get caught while they are still cheap to fix.

### Alternatives considered
- **TimeTree.org as the primary tree source instead of the VertLife/Upham mammal supertree.** TimeTree is easier to query and spans the whole tree of life, not just mammals, which makes it a reasonable default for quick lookups. It wasn't chosen as the primary source here because it typically returns a single point-estimate tree, whereas Upham et al.'s mammal-specific supertree ships as a distribution of trees capturing topological and branch-length uncertainty — this lets Phase 5's PGLS models be re-run across multiple trees from that distribution to check whether the results are sensitive to phylogenetic uncertainty, a check TimeTree's single-tree output doesn't support. For an all-mammal study, Upham et al. also has denser within-mammal taxon sampling than TimeTree's more general-purpose coverage.
- **Build a new tree from scratch by aligning marker genes across the 55 overlap species.** This would give full control over the exact species set and calibration choices. It wasn't chosen because it's a substantial standalone bioinformatics undertaking (ortholog selection, multiple sequence alignment, tree search, divergence-time calibration) that essentially re-derives something already solved: existing, peer-reviewed, actively maintained trees (TimeTree, VertLife) already cover essentially all 55 overlap species, so building a new tree would spend significant effort duplicating existing, validated work rather than advancing the actual hypothesis test.

### Deliverable
A pruned, branch-length-resolved phylogeny file (Newick or NEXUS) containing exactly the Phase 1 overlap species, plus a note of any species dropped at this stage and why.

---

## Phase 3 — Bioinformatics processing of the microbiome data

**Objective:** Turn raw or partially processed sequencing data into quality-controlled, taxonomically classified, genome-resolved microbial data (MAGs).

### Steps

1. **Quality control and preprocessing.** Trim sequencing adapters and low-quality bases, then remove host DNA reads. `fastp` is the current standard single-pass tool for adapter/quality trimming — it combines what used to require separate tools (FastQC-style QC reporting plus Trimmomatic-style trimming) into one fast pass over the data [12]. Host-read removal requires mapping reads against each sample's *own host species* reference genome (or nearest available relative) and discarding matches — a nonstandard wrinkle for a 55-species, multi-order comparative dataset like this one, versus the single-reference-genome case (e.g., always mapping against the human genome) that most metagenomics pipeline documentation assumes.

2. **Quick taxonomic profiling before committing to full assembly.** Run a k-mer-based classifier (Kraken2 [13]) or marker-gene profiler (MetaPhlAn 4 [14]) on the QC'd reads as an early sanity check on community composition — this is much faster than assembly and catches sample-swap or contamination problems before expensive compute is spent on a bad sample.

3. **Assembly and binning into MAGs.** This is the expensive step, and the recommended path is a maintained, containerized pipeline rather than hand-chaining individual tools: `nf-core/mag` is a Nextflow-based, community-maintained pipeline that runs QC, assembly, and binning end-to-end and is designed to be reproducible across compute environments [15]. Quality-filter the resulting bins with CheckM2, a machine-learning-based genome completeness/contamination estimator that has replaced the original CheckM as the current standard [16], and assign taxonomy with GTDB-Tk, which classifies genomes against the Genome Taxonomy Database rather than older, less internally consistent naming schemes [17].

4. **Reuse published MAGs where they already exist, rather than re-assembling.** Youngblut et al. 2020 already published 5,596 assembled genomes collapsing to 1,522 species-level genome bins (SGBs) from 180 host species processed with this same general workflow [4] — for any Phase 1 overlap species covered by that dataset, this step is already done, which is a strong practical argument for prioritizing species covered by Youngblut et al. when the overlap allows a choice.

### Best-practice rationale
Bharti and Grimm's 2021 review of microbiome analysis best practices is the most current, tool-agnostic reference for why each of these stages exists and in what order: preprocessing/QC must precede any downstream step because contamination and adapter content otherwise inflate false taxonomic diversity and corrupt assembly graphs; genome-quality filtering (CheckM2) must happen before any genome is used for functional inference, because a 40%-complete "genome" produces unreliable gene-content calls [11]. Using a single maintained pipeline (nf-core/mag) rather than hand-assembled scripts is recommended specifically for reproducibility — the pipeline documents its own software versions and parameters as part of its output, which matters when the project is later revisited or reviewed [15].

### Alternatives considered
- **QC/trimming: `fastp` vs. Trimmomatic [38].** Trimmomatic was the long-standing standard and is highly configurable via explicit rule chains, but published benchmarks show `fastp` runs several-fold faster with comparable trimming quality, and it generates its own QC report in the same pass, removing the need for a separate FastQC step. Given the number of samples this project needs to process (up to ~103 mammal-species samples if the full Youngblut cohort is used), the throughput advantage is the deciding factor.
- **Taxonomic sanity check: Kraken2 vs. MetaPhlAn4.** Kraken2's k-mer-based approach classifies a larger fraction of reads (higher sensitivity, and can go to strain level when paired with Bracken) but is more prone to false-positive calls driven by reference-database contamination; MetaPhlAn4 uses a curated marker-gene set for higher-precision but lower-recall species calls. Because this step in Phase 3 is only a fast early sanity check — not the dataset's final abundance table, which comes from the MAG-based work in Phase 4 — the higher-sensitivity Kraken2 is preferred here specifically for quickly catching swapped or contaminated samples; a MetaPhlAn4 pass is worth revisiting later if marker-gene-based relative abundance is wanted as an independent cross-check on the genome-based feature calls.
- **Assembly/binning: `nf-core/mag` vs. a hand-chained pipeline (e.g., SPAdes + MetaBAT2 [39]) or the anvi'o interactive platform [40].** anvi'o is well suited to hands-on, interactive bin refinement and visualization, which is valuable when a person wants to manually inspect and curate each bin — a reasonable choice for a small, high-value sample set. It wasn't chosen here because this project needs a reproducible, versioned pipeline that scales across potentially 100+ samples without manual per-sample intervention; `nf-core/mag`'s Nextflow structure automatically logs the exact tool versions and parameters used, which satisfies the Phase 8 reproducibility goal without separate manual bookkeeping.
- **Genome quality filtering: CheckM2 vs. the original CheckM.** CheckM1 estimates completeness/contamination using lineage-specific single-copy marker gene sets, which assumes the genome's lineage is already well represented in reference databases. Many gut MAGs recovered from non-model host species (most of the animals in this project) belong to novel, poorly characterized microbial lineages, which is exactly where CheckM1's lineage-specific approach is known to degrade. CheckM2's machine-learning-based approach was trained and benchmarked to be more accurate on such novel/undersampled lineages, which is why it's used instead.
- **Taxonomy assignment: GTDB-Tk vs. NCBI taxonomy.** NCBI taxonomy is more familiar and human-readable, but it is known to retain polyphyletic and inconsistently ranked groups because it partly reflects historical naming conventions rather than current phylogenetic evidence. GTDB-Tk assigns taxonomy from genome-based phylogenetic placement, producing internally consistent species-level clusters — which matters here specifically for deciding whether MAGs recovered from two different host species represent the "same" or "different" microbial species when comparing functional content across the dataset.

### Deliverable
Per-sample QC reports, a taxonomic profile table, and a set of quality-filtered, taxonomically classified MAGs for each Phase 1 overlap species (reused from Youngblut et al. where possible, newly generated otherwise).

---

## Phase 4 — Functional annotation and feature matrix construction

**Objective:** Convert genome sequence into a per-species table of functional features that can serve as the independent-variable matrix for Phase 5.

### Steps

1. **Annotate each MAG's genes against a metabolic pathway database.** KEGG Orthology (KO) assignment via a profile-HMM tool (KofamScan/KofamKOALA [19]) or eggNOG-mapper (which additionally gives broader orthology and COG functional categories in the same pass) [18] are the two most common current choices; eggNOG-mapper is particularly suited to metagenome-scale annotation because it was explicitly built and benchmarked for that use case [18].

2. **Annotate carbohydrate-active enzyme content separately.** CAZymes (glycoside hydrolases, glycosyltransferases, etc.) are not well captured by general KO annotation and need a dedicated database; dbCAN3 is the current maintained tool, and as of its 2023 version adds substrate-level prediction on top of enzyme-family assignment, which is directly useful for asking "what could this organism actually break down" rather than just "does it carry a carbohydrate-active enzyme" [20].

3. **Collapse gene-level annotations into a per-species feature matrix, not a per-gene or per-MAG one.** The unit of analysis in Phase 5 is the species (to match the cancer-mortality data), so gene/pathway presence needs to be aggregated across all MAGs assigned to a given host species — e.g., pathway presence/absence, gene family relative abundance, or specific gene counts (bile-acid-modifying genes, butyrate-synthesis genes, colibactin-associated *pks* genes) per species.

4. **Do not treat this feature matrix as ordinary continuous data without transformation.** Microbiome-derived abundance and count data are compositional — each sample sums to an arbitrary total set by the sequencing instrument, not by the underlying biology — and applying standard statistical methods (correlation, standard regression) to raw proportions or counts produces spurious results specifically because of this constraint. The standard fix is a log-ratio transformation (most commonly centered log-ratio, CLR) before any downstream statistical modeling [21]. This step is easy to skip because the data "looks like" ordinary continuous data; it is not.

### Best-practice rationale
Annotating against two purpose-built databases (KEGG/eggNOG for general pathways, dbCAN for carbohydrate metabolism specifically) rather than one general-purpose annotation captures functional categories that matter directly for this project's downstream questions (fiber fermentation → SCFA production is a CAZyme-driven process, not a generic KO-driven one) [18,20]. The compositionality correction is not optional despite frequently being skipped in practice — the paper that established this as a hard requirement rather than a stylistic preference is titled, deliberately, "Microbiome Datasets Are Compositional: And This Is Not Optional" [21].

### Alternatives considered
- **Annotate with only eggNOG-mapper, or only KofamKOALA, instead of both.** Running a single tool is faster and simpler to maintain. Both are used here instead because they rely on different underlying HMM profile sets and can disagree at the margins for any given gene; cross-validating hits between the two specifically for the biologically important gene set this project cares about most (bile-acid, butyrate, and colibactin pathway genes) catches annotation-tool-specific false positives for exactly the genes the downstream mechanistic story depends on — worth the added runtime for that subset even though it isn't run genome-wide for every gene.
- **CAZyme annotation via a manual BLAST search against the CAZy database, instead of dbCAN3.** Manual BLAST against CAZy was the older standard approach, but it requires hand-curated e-value/identity thresholds per enzyme family and gives no substrate-level prediction. dbCAN3 combines three independent detection methods (HMMER, DIAMOND, and dbCAN_sub) into one consensus call and adds substrate-level prediction on top, giving both higher-confidence family assignments and the "what can this actually break down" information this project's fiber-fermentation question specifically needs.
- **Rarefaction (subsampling all samples to equal read depth) instead of CLR transformation.** Rarefaction was historically the standard fix for uneven sequencing depth across samples. It's not used here because it has been shown to discard a large fraction of valid read data without actually resolving the underlying compositionality problem, and is now explicitly recommended against for this kind of feature-comparison analysis [41]. CLR transformation addresses compositionality directly, without throwing away data, which is why it's used instead.

### Deliverable
A species × feature matrix (KEGG pathway presence/abundance, CAZy family abundance, and specific gene-of-interest counts), CLR-transformed and ready to serve as the independent variable set in Phase 5.

---

## Phase 5 — Statistical correlation with phylogenetic correction

**Objective:** Test whether microbiome features (Phase 4) predict cancer mortality (Phase 1), correcting for the fact that species are not statistically independent data points.

### Steps

1. **Fit PGLS, not ordinary regression.** Phylogenetic Generalized Least Squares replaces the identity error-covariance matrix assumed by ordinary least squares with one derived from the Phase 2 tree, so that closely related species (who share both traits and ancestry) don't produce inflated, spurious correlations [23,24]. In R, this is implemented in the `caper` package (`pgls()` function) [25] or the faster `phylolm` package, which is more practical when testing many features (e.g., hundreds of KEGG pathways) because its underlying algorithm scales linearly rather than cubically with the number of species [28].

2. **Estimate the phylogenetic signal parameter (Pagel's λ) rather than fixing it.** λ ranges from 0 (no phylogenetic structure — equivalent to ordinary regression) to 1 (full Brownian-motion phylogenetic structure); both `caper` and `phylolm` can estimate λ by maximum likelihood alongside the regression coefficients rather than requiring the analyst to assume a value, and doing so is the standard recommended practice — assuming λ=1 when the real data show weaker phylogenetic signal, or vice versa, biases the resulting standard errors [24,26].

3. **Include the covariates fixed in Phase 0, in the same model, not as separate univariate tests.** Litter size (from Boddy et al. [2]) and diet category (from Milani et al. [3]) are known or plausible confounders of both microbiome composition and cancer risk; a univariate PGLS of cancer mortality against a single microbiome feature, run without these covariates, cannot distinguish a real microbiome effect from a shared-confound artifact.

4. **Correct for multiple testing across all features tested in the same family.** If testing hundreds of KEGG pathways or gene families, apply the Benjamini-Hochberg false discovery rate procedure across that full set of tests rather than reading off raw p-values — this is the standard correction for exactly this "many tests off one dataset" scenario and controls the *expected proportion* of false positives among the features called significant, rather than the stricter (and often unnecessarily conservative for exploratory work) family-wise error rate [22].

5. **Check the residuals and flag phylogenetic outliers.** Species whose observed cancer mortality is very different from what the fitted PGLS model (given their microbiome features and phylogenetic position) predicts are worth flagging individually — they may represent genuinely interesting biology (a species that has evolved unusual cancer resistance despite an otherwise "risky" microbiome profile) or a data-quality problem, and either way are worth a manual look before being folded silently into a p-value.

### Best-practice rationale
Ignoring phylogenetic non-independence in comparative datasets was identified decades ago as producing an inflated, unreliable rate of false-positive associations, which is the entire motivation for PGLS existing as a method [23,24]. The Symonds & Blomberg primer [26] and Mundry's companion chapter on PGLS's statistical assumptions [27] are the standard practical starting references for actually running this correctly rather than just conceptually — Mundry's chapter in particular focuses on exactly the kind of small-sample, real-data pitfalls (non-normality, heteroscedasticity, outliers) this project's ~55-species dataset is likely to hit [27]. This step is also the direct statistical analog of the Youngblut et al. 2019 diet/phylogeny study discussed earlier in this project [29], which used PGLS to separate diet's effect on microbiome diversity from phylogeny's effect using the same logic applied here to cancer mortality.

### Alternatives considered
- **Felsenstein's phylogenetic independent contrasts (PIC) [42] instead of PGLS.** PIC was the original solution to phylogenetic non-independence and is the historical precursor to PGLS. It wasn't chosen because it only handles a single continuous predictor cleanly under a strict Brownian-motion model of evolution, and doesn't readily accommodate multiple covariates in one model. PGLS is a generalization of the same underlying idea within a standard generalized-least-squares regression framework, which supports the multiple covariates (litter size, diet) Phase 0 already fixed, plus flexible correlation structures via the estimated λ — which is why PGLS, not PIC, is used here.
- **Bayesian phylogenetic mixed models (e.g., MCMCglmm) [43] instead of frequentist PGLS.** MCMCglmm-style models can jointly integrate over both phylogenetic uncertainty and tree topology uncertainty (useful given Phase 2's multi-tree VertLife output) and are often considered the more statistically thorough option. They weren't chosen as the primary screening tool because they require substantially more setup per model (prior specification, MCMC convergence diagnostics) and computation time, which isn't practical when screening potentially hundreds of candidate pathway features. The frequentist PGLS/`phylolm` approach is used for that large-scale screen; a Bayesian mixed-model re-analysis is a reasonable confirmatory follow-up for the short list of features that survive the Phase 5 FDR correction, rather than the tool used for the initial pass.
- **`caper` vs. `phylolm` as the specific PGLS implementation.** Both are used, for different purposes: `caper`'s `pgls()` gives more complete model diagnostics (residual plots, profile likelihood for λ) and is used for the final, small number of confirmatory models on FDR-surviving features; `phylolm` is used for the large-scale feature screen because its underlying algorithm scales linearly, not cubically, with the number of species, which matters when fitting potentially hundreds of models across a tree of this size — using `caper` for that many models would be needlessly slow without added benefit at the screening stage.
- **Bonferroni correction or permutation-based correction instead of Benjamini-Hochberg FDR.** Bonferroni controls the family-wise error rate very conservatively; given the project may test hundreds of correlated pathway features (many pathways share genes and co-vary), Bonferroni would likely eliminate almost all true positives along with the false ones. Permutation-based correction (shuffling trait labels across the tree) handles feature correlation structure better but is far more computationally expensive to run per feature at the exploratory-screening stage used here. Benjamini-Hochberg FDR is the standard middle ground for this kind of exploratory, many-features-one-dataset scenario, with permutation-based correction left as an option for the final confirmatory model on the short list of FDR-surviving hits.

### Deliverable
A table of PGLS results (effect size, t-value, p-value, FDR-adjusted p-value, estimated λ) for each tested microbiome feature, with covariates included, plus a flagged list of phylogenetic outlier species.

---

## Phase 6 — Mechanistic follow-up on the hits

**Objective:** For features that survive Phase 5's statistical test, build a specific, checkable mechanistic story rather than stopping at "this pathway correlates with cancer mortality."

### Steps

1. **Build genome-scale metabolic models (GEMs) from the relevant MAGs.** CarveMe is the standard tool for fast, automated GEM reconstruction directly from a genome/MAG, including a mode built specifically for constructing community-level (multi-species) models rather than only single-organism ones [30]. `gapseq` is a commonly used alternative/complement, notably strong at pathway-level prediction and gap-filling — reconstructing plausible reactions to complete pathways that are only partially present in the genome annotation, which matters because MAG annotations are frequently incomplete [31].

2. **Run flux balance analysis (FBA) under a diet-informed nutrient constraint, not a generic/default medium.** FBA predicts the feasible range of metabolic flux through a genome-scale network given a set of nutrient inputs; COBRApy is the standard Python toolkit for actually running this optimization once a model exists [32]. Using the actual diet category recorded for each host species (already present in the Phase 1 microbiome datasets, e.g., herbivore/carnivore/omnivore) as the nutrient-input constraint, rather than a generic laboratory growth medium, is what makes the resulting flux predictions biologically meaningful for this specific comparative question.

3. **Model the community, not just isolated single species.** Real gut communities cross-feed — one organism's metabolic byproduct is another's substrate — so single-organism FBA will miss interaction-dependent metabolite production. MICOM is a purpose-built framework for exactly this: metagenome-scale community modeling that has been validated against observed *in vivo* growth rates in gut communities [33].

4. **Cross-reference predicted output against a curated list of cancer-relevant metabolites.** Secondary bile acids (produced by microbial 7α-dehydroxylation of primary bile acids), hydrogen sulfide, and genotoxins like colibactin are established or strongly suspected pro-carcinogenic bacterial metabolites; short-chain fatty acids (particularly butyrate) are established protective ones, largely through anti-inflammatory and anti-proliferative effects on colonic epithelium — this is reviewed comprehensively in Louis, Hold & Flint's 2014 Nature Reviews Microbiology synthesis of gut-microbiota-derived metabolites and colorectal cancer [34]. Checking whether the Phase 6 FBA output for a Phase-5-significant species includes elevated predicted flux through any of these specific compounds turns a statistical association into a testable, literature-grounded mechanistic hypothesis.

### Best-practice rationale
FBA predictions are only as trustworthy as the constraints they're run under — an FBA model run on a default/generic medium will produce a "what this organism could theoretically do under lab conditions" answer, not "what this organism is actually doing in this animal's gut," which is why diet-constrained, community-level modeling (rather than single-species, default-medium modeling) is emphasized here specifically [30,33]. This phase is explicitly framed as *mechanistic follow-up*, not primary hypothesis testing — its job is to generate a specific, falsifiable claim ("species X's community is predicted to overproduce deoxycholic acid") that Phase 7 can then actually test.

### Alternatives considered
- **`gapseq` as the primary bulk GEM-reconstruction tool instead of CarveMe.** `gapseq`'s pathway-level, database-informed gap-filling tends to produce more biologically accurate models for specific pathways — including the anaerobic fermentation and bile-acid pathways this project cares most about — but it's markedly slower and less practical for reconstructing models in bulk across potentially dozens of MAGs per species. CarveMe is used as the primary bulk-reconstruction tool for its speed and dedicated community-model mode, with `gapseq` reserved for a targeted second-pass rebuild of the specific bile-acid/butyrate/colibactin pathway regions in the Phase-5 hit species, where the extra per-pathway accuracy is worth the added runtime.
- **SteadyCom [44] instead of MICOM for community-level flux modeling.** SteadyCom is a well-established alternative that also predicts steady-state community composition and fluxes. It wasn't chosen because its reference implementation runs through the MATLAB-based COBRA Toolbox, which doesn't fit this project's existing Python/COBRApy tooling. MICOM is Python-native, integrates directly with COBRApy models produced by CarveMe/`gapseq`, and has been validated specifically against observed in vivo gut community growth rates — all of which fit this project's stack and validation needs better.
- **Simple gene-presence/absence scoring instead of full FBA.** Counting bile-acid-modifying, butyrate-synthesis, or *pks* (colibactin) genes per species directly from the Phase 4 annotation table is much cheaper and was in fact already partly explored as part of Phase 4. It's not sufficient on its own here because gene presence doesn't guarantee pathway flux — a gene can be present but flux-limited by missing cofactor-supply reactions elsewhere in the network, or the reaction may not be able to carry flux under the organism's actual nutrient constraints. FBA is used in Phase 6 specifically to move from "has the gene" to "the network is predicted to actually be capable of the flux," a stronger and more specific mechanistic claim than presence/absence alone can support.

### Deliverable
Genome-scale community metabolic models for the Phase 5 hit species, FBA-predicted metabolite output under each species' actual diet constraint, and a short list of specific compounds flagged for downstream validation.

---

## Phase 7 — Experimental validation (longer-term, separate proposal)

**Objective:** Move from correlation/prediction toward causal evidence for the strongest Phase 6 candidates. This phase is realistically a follow-on study, not part of the initial computational analysis, but should be scoped now so the proposal states an endpoint rather than stopping at "we found a correlation."

### Steps

1. **Select candidate taxa/pathways, not the full hit list, for validation.** Experimental validation is slow and expensive relative to everything upstream; prioritize the Phase 6 candidates with the clearest mechanistic story (a specific predicted metabolite with an established cancer-relevant mechanism) over marginal statistical hits.

2. **Use fecal microbiota transplant (FMT) into germ-free (gnotobiotic) mice as the standard causal-inference design.** Colonizing a microbiome-naive host with a donor community and observing the resulting phenotype is the standard way comparative microbiome findings get tested for causality, comparing transplants from high-cancer-mortality-associated versus low-cancer-mortality-associated source microbiomes.

3. **Design the experiment with the field's known failure modes for this exact design in mind, not just standard animal-study rigor.** A systematic review of published human-microbiota-associated (HMA) rodent studies found that 95% reported successful transfer of the donor phenotype to the recipient animal — a rate the authors argue is implausibly high and likely reflects publication bias, small effect sizes being over-interpreted, and confounding from incomplete or unstable colonization, rather than the gut microbiome being that reliably causal across this many independent studies [35]. Concretely, this means: pre-register the specific phenotype being tested (extending Phase 0's practice to this experiment), use adequate sample sizes per group (informed by that paper's discussion of typical underpowering in this literature), and report negative or partial-transfer results rather than only positive ones.

4. **Report the experiment using the ARRIVE 2.0 checklist.** This is the current standard structure for reporting in vivo animal research — covering study design, sample size justification, randomization, and outcome reporting — and using it from the experimental design stage (not just retrofitted at write-up) makes the eventual results easier to evaluate and reproduce [36].

### Best-practice rationale
Gnotobiotic FMT is the accepted causal-inference tool in this field, but the Walter et al. 2020 Cell paper is included here specifically as a corrective: it is the most direct, current source on why this exact experimental design over-claims causality more often than it should, and what to do differently [35]. Citing it here is meant to set expectations correctly in the proposal — that a positive FMT result would be suggestive, not automatically definitive, and that this needs to be stated in the write-up rather than discovered by a reviewer.

### Alternatives considered
- **In vitro gut-simulator systems (e.g., SHIME-type multi-stage bioreactors) instead of, or before, germ-free mouse FMT.** In vitro bioreactor systems are cheaper, faster, and avoid animal use, and are a reasonable early screen for whether a donor community actually produces the Phase 6-predicted metabolites in a simplified gut-like environment. They aren't sufficient as the final validation step because they can't capture host-immune or host-epithelial interactions, which is where the cancer-relevant phenotype (tumorigenesis, inflammation) actually plays out. A bioreactor pre-screen is a reasonable, cheaper addition before committing to the mouse study, but germ-free FMT remains the necessary endpoint for an actual causal cancer-phenotype claim.
- **Antibiotic-depletion mouse models instead of germ-free (gnotobiotic) mice.** Depleting an existing mouse gut community with broad-spectrum antibiotics before introducing the donor community is cheaper and doesn't require specialized gnotobiotic animal facilities. It wasn't chosen because it leaves a residual, incompletely eliminated community and residual antibiotic effects on host physiology, both of which confound interpretation of the transplanted phenotype. Germ-free mice avoid this confound entirely, which is why they're specified here despite the higher facility cost — consistent with Walter et al.'s [35] critique of exactly this kind of incomplete-colonization confound as a driver of over-claimed causality in the HMA-rodent literature.

### Deliverable
A short, separately scoped experimental proposal: candidate taxa/metabolites selected from Phase 6, an FMT + germ-free mouse design with pre-specified phenotype and sample size, and a reporting plan following ARRIVE 2.0 [36].

---

## Phase 8 — Write-up

**Objective:** Report the findings in a way that states their real strength and limitations up front, rather than requiring a reviewer to discover them.

### Steps

1. **Report effect sizes and confidence intervals for every tested feature, not just p-values.** A p-value alone conveys whether an effect was detected at some threshold, not how large or biologically meaningful it is; standardized effect-size statistics (e.g., correlation-type *r* or standardized regression coefficients) with confidence intervals let a reader judge both the magnitude and the precision of each estimate, and are the recommended standard reporting format for exactly this reason [37].

2. **State the sample-size/power limitation explicitly, before results are presented, not after.** With a species-level *n* in the tens (55 species in this project's current overlap) tested against potentially hundreds of pathway-level features, this design has drastically less statistical power than the GWAS-style designs it is methodologically analogous to (see the earlier GWAS-comparison discussion in this project) — GWAS studies achieve power with sample sizes in the thousands-to-millions specifically because per-variant effect sizes are small and the multiple-testing burden is severe; this project's much smaller *n* means it can realistically detect only comparatively large effects, and many true but modest associations will be missed (Type II error), not just falsely detected ones controlled for by the Phase 5 FDR correction. Say this in the introduction or methods, not the discussion's limitations paragraph as an afterthought.

3. **Report the phylogenetic effective sample size alongside the raw species count.** Because closely related species are not independent data points, the raw count of species overstates the amount of independent information the dataset actually contains; reporting the estimated λ from Phase 5 alongside the species count gives a reader a more honest picture of how much phylogenetically independent signal the dataset really has.

4. **Distinguish, explicitly, which findings are Phase 5 statistical associations versus Phase 6 mechanistic predictions versus (if reached) Phase 7 experimental results**, rather than letting the strength of causal language drift upward between sections. A PGLS association and an FBA-predicted metabolite are both interesting, but neither is causal evidence on its own — only Phase 7 addresses causality directly, and even then with the caveats noted above [35].

### Best-practice rationale
Nakagawa and Cuthill's practical guide is the standard reference specifically for why effect sizes and confidence intervals, not p-values alone, are the appropriate unit of biological reporting — null-hypothesis significance testing tells a reader whether an effect was detectable, not how large or important it is [37]. Stating power/sample-size limitations upfront rather than as a defensive afterthought is standard practice in any comparative study with a small, hard-capped *n* (here, capped by how many species happen to have both cancer-mortality and microbiome data available, not by study design choice) — and is exactly the kind of limitation this project's GWAS comparison earlier in the process was meant to surface early.

### Alternatives considered
- **Report only p-values with a significance threshold (the more common historical default in comparative biology write-ups).** This requires less additional computation and is still what many readers expect. It isn't used as the primary reporting format here because, as covered by Nakagawa & Cuthill [37], p-values alone conflate statistical detectability with biological importance and don't communicate how precise or imprecise an estimate is — which matters more than usual for a project with this small an *n*. Effect sizes with confidence intervals are reported instead specifically to make that uncertainty visible rather than collapsing it into a binary "significant / not significant" call.

### Deliverable
A written report/manuscript draft with: hypothesis and covariates (from Phase 0), species overlap and data provenance (Phase 1), tree and pruning notes (Phase 2), bioinformatics pipeline and QC summary (Phase 3), feature matrix construction (Phase 4), full PGLS results table with effect sizes and CIs (Phase 5), FBA-based mechanistic candidates (Phase 6), experimental proposal if applicable (Phase 7), and an explicit power/sample-size limitations statement placed in the methods, not buried in the discussion.

---

## What else would help move this forward

A few concrete things that would sharpen or accelerate specific phases, if available:

- **Raw sequencing data (FASTQ) or SRA/ENA accession numbers** for any of the 55 overlap species not already covered by Youngblut et al.'s published MAGs [4] — this is the single biggest remaining gap between "data acquired" (Phase 1, done) and "bioinformatics processing" (Phase 3, not yet started).
- **A decision on the Phase 0 hypothesis statement** (taxonomic vs. functional vs. metabolite-level endpoint) — the rest of this document assumes the metabolite-level endpoint per the project's stated FBA interest, but this should be confirmed rather than assumed.
- **Any existing diet-composition data** (specific nutrient/macronutrient breakdowns, not just diet category) for the overlap species, which would sharpen the Phase 6 FBA diet constraints beyond the coarse herbivore/carnivore/omnivore categories currently available.
- **Access to a compute environment with full network access** (this sandbox's network is allow-listed and currently blocks NCBI, EBI, MassIVE, Zenodo, and Dryad, which is why several Phase 1 sources needed to be manually downloaded rather than auto-fetched) — needed for Phase 3's SRA/ENA pulls and any GBIF/NCBI Taxonomy synonym-resolution work in Phase 1/2.

---

## References

[1] Vincze O, Colchero F, Lemaitre JF, Conde DA, Pavard S, Bieuville M, Urrutia AO, Ujvari B, Boddy AM, Maley CC, Thomas F, Giraudeau M. 2022. Cancer risk across mammals. *Nature* 601:263–267. https://doi.org/10.1038/s41586-021-04224-5

[2] Boddy AM, Abegglen LM, Pessier AP, Aktipis A, Schiffman JD, Maley CC, Witte C. 2020. Lifetime cancer prevalence and life history traits in mammals. *Evolution, Medicine, and Public Health* 2020(1):187–195. https://doi.org/10.1093/emph/eoaa015

[3] Milani C, Alessandri G, Mancabelli L, Mangifesta M, Lugli GA, Viappiani A, Longhi G, Anzalone R, Duranti S, Turroni F, Ossiprandi MC, van Sinderen D, Ventura M. 2020. Multi-omics Approaches To Decipher the Impact of Diet and Host Physiology on the Mammalian Gut Microbiome. *Applied and Environmental Microbiology* 86:e01864-20. https://doi.org/10.1128/AEM.01864-20

[4] Youngblut ND, de la Cuesta-Zuluaga J, Reischer GH, Dauser S, Schuster N, Walzer C, Stalder G, Farnleitner AH, Ley RE. 2020. Large-Scale Metagenome Assembly Reveals Novel Animal-Associated Microbial Genomes, Biosynthetic Gene Clusters, and Other Genetic Diversity. *mSystems* 5:e01045-20. https://doi.org/10.1128/mSystems.01045-20

[5] Gregor R, Probst M, Eyal S, Aksenov A, Sasson G, Horovitz I, Dorrestein PC, Meijler MM, Mizrahi I. 2022. Mammalian gut metabolomes mirror microbiome composition and host phylogeny. *The ISME Journal* 16:1262–1274. https://doi.org/10.1038/s41396-021-01152-0

[6] Nosek BA, Ebersole CR, DeHaven AC, Mellor DT. 2018. The preregistration revolution. *Proceedings of the National Academy of Sciences* 115(11):2600–2606. https://doi.org/10.1073/pnas.1708274114

[7] Chamberlain SA, Szöcs E. 2013. taxize: taxonomic search and retrieval in R. *F1000Research* 2:191. https://doi.org/10.12688/f1000research.2-191.v2

[8] Kumar S, Suleski M, Craig JM, Kasprowicz AE, Sanderford M, Li M, Stecher G, Hedges SB. 2022. TimeTree 5: An Expanded Resource for Species Divergence Times. *Molecular Biology and Evolution* 39(8):msac174. https://doi.org/10.1093/molbev/msac174

[9] Upham NS, Esselstyn JA, Jetz W. 2019. Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. *PLoS Biology* 17(12):e3000494. https://doi.org/10.1371/journal.pbio.3000494

[10] Bharti R, Grimm DG. 2021. Current challenges and best-practice protocols for microbiome analysis. *Briefings in Bioinformatics* 22(1):178–193. https://doi.org/10.1093/bib/bbz155

[11] (see [10])

[12] Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884–i890. https://doi.org/10.1093/bioinformatics/bty560

[13] Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. *Genome Biology* 20:257. https://doi.org/10.1186/s13059-019-1891-0

[14] Blanco-Míguez A, Beghini F, Cumbo F, McIver LJ, Thompson KN, Zolfo M, et al. 2023. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology* 41(11):1633–1644. https://doi.org/10.1038/s41587-023-01688-w

[15] Krakau S, Straub D, Gourlé H, Gabernet G, Nahnsen S. 2022. nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning. *NAR Genomics and Bioinformatics* 4(1):lqac007. https://doi.org/10.1093/nargab/lqac007

[16] Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods* 20:1203–1212. https://doi.org/10.1038/s41592-023-01940-w

[17] Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. 2022. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics* 38(23):5315–5316. https://doi.org/10.1093/bioinformatics/btac672

[18] Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. 2021. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. *Molecular Biology and Evolution* 38(12):5825–5829. https://doi.org/10.1093/molbev/msab293

[19] Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, Ogata H. 2020. KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics* 36(7):2251–2252. https://doi.org/10.1093/bioinformatics/btz859

[20] Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y. 2023. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research* 51(W1):W115–W121. https://doi.org/10.1093/nar/gkad328

[21] Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. 2017. Microbiome Datasets Are Compositional: And This Is Not Optional. *Frontiers in Microbiology* 8:2224. https://doi.org/10.3389/fmicb.2017.02224

[22] Benjamini Y, Hochberg Y. 1995. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society Series B* 57(1):289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

[23] Pagel M. 1999. Inferring the historical patterns of biological evolution. *Nature* 401:877–884. https://doi.org/10.1038/44766

[24] Freckleton RP, Harvey PH, Pagel M. 2002. Phylogenetic Analysis and Comparative Data: A Test and Review of Evidence. *The American Naturalist* 160(6):712–726. https://doi.org/10.1086/343873

[25] Orme D, Freckleton R, Thomas G, Petzoldt T, Fritz S, Isaac N, Pearse W. 2013. caper: Comparative Analyses of Phylogenetics and Evolution in R. R package version 0.5.2. https://cran.r-project.org/package=caper

[26] Symonds MRE, Blomberg SP. 2014. A primer on phylogenetic generalised least squares. In: Garamszegi LZ (ed) *Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology*, Chapter 5, pp. 105–130. Springer, Berlin. https://doi.org/10.1007/978-3-662-43550-2_5

[27] Mundry R. 2014. Statistical Issues and Assumptions of Phylogenetic Generalized Least Squares. In: Garamszegi LZ (ed) *Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology*, Chapter 6. Springer, Berlin. https://doi.org/10.1007/978-3-662-43550-2_6

[28] Ho LST, Ané C. 2014. A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Systematic Biology* 63(3):397–408. https://doi.org/10.1093/sysbio/syu005

[29] Youngblut ND, Reischer GH, Walters W, Schuster N, Walzer C, Stalder G, Ley RE, Farnleitner AH. 2019. Host diet and evolutionary history explain different aspects of gut microbiome diversity among vertebrate clades. *Nature Communications* 10:2200. https://doi.org/10.1038/s41467-019-10191-3

[30] Machado D, Andrejev S, Tramontano M, Patil KR. 2018. Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Research* 46(15):7542–7553. https://doi.org/10.1093/nar/gky537

[31] Zimmermann J, Kaleta C, Waschina S. 2021. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. *Genome Biology* 22:81. https://doi.org/10.1186/s13059-021-02295-1

[32] Ebrahim A, Lerman JA, Palsson BO, Hyduke DR. 2013. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Systems Biology* 7:74. https://doi.org/10.1186/1752-0509-7-74

[33] Diener C, Gibbons SM, Resendis-Antonio O. 2020. MICOM: Metagenome-Scale Modeling To Infer Metabolic Interactions in the Gut Microbiota. *mSystems* 5(1):e00606-19. https://doi.org/10.1128/mSystems.00606-19

[34] Louis P, Hold GL, Flint HJ. 2014. The gut microbiota, bacterial metabolites and colorectal cancer. *Nature Reviews Microbiology* 12:661–672. https://doi.org/10.1038/nrmicro3344

[35] Walter J, Armet AM, Finlay BB, Shanahan F. 2020. Establishing or Exaggerating Causality for the Gut Microbiome: Lessons from Human Microbiota-Associated Rodents. *Cell* 180(2):221–232. https://doi.org/10.1016/j.cell.2019.12.025

[36] Percie du Sert N, Hurst V, Ahluwalia A, Alam S, Avey MT, Baker M, et al. 2020. The ARRIVE guidelines 2.0: Updated guidelines for reporting animal research. *PLoS Biology* 18(7):e3000410. https://doi.org/10.1371/journal.pbio.3000410

[37] Nakagawa S, Cuthill IC. 2007. Effect size, confidence interval and statistical significance: a practical guide for biologists. *Biological Reviews* 82(4):591–605. https://doi.org/10.1111/j.1469-185X.2007.00027.x

[38] Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* 30(15):2114–2120. https://doi.org/10.1093/bioinformatics/btu170

[39] Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. 2019. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7:e7359. https://doi.org/10.7717/peerj.7359

[40] Eren AM, Esen ÖC, Quince C, Vineis JH, Morrison HG, Sogin ML, Delmont TO. 2015. Anvi'o: an advanced analysis and visualization platform for 'omics data. *PeerJ* 3:e1319. https://doi.org/10.7717/peerj.1319

[41] McMurdie PJ, Holmes S. 2014. Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible. *PLoS Computational Biology* 10(4):e1003531. https://doi.org/10.1371/journal.pcbi.1003531

[42] Felsenstein J. 1985. Phylogenies and the Comparative Method. *The American Naturalist* 125(1):1–15. https://doi.org/10.1086/284325

[43] Hadfield JD. 2010. MCMC Methods for Multi-Response Generalized Linear Mixed Models: The MCMCglmm R Package. *Journal of Statistical Software* 33(2):1–22. https://doi.org/10.18637/jss.v033.i02

[44] Chan SHJ, Simons MN, Maranas CD. 2017. SteadyCom: Predicting microbial abundances while ensuring community stability. *PLoS Computational Biology* 13(5):e1005539. https://doi.org/10.1371/journal.pcbi.1005539

---

*A note on data provenance in this repository:* the microbiome dataset referred to internally in `phase1_data_ingestion.py` as `manor2020_microbiome` is authored by Milani et al. [3], not "Manor" — that variable name was an incorrect guess made before the actual author list was verified against the published paper. It's a naming issue only (the underlying data and citation in this document are correct); renaming the variable/file is a small cleanup worth doing before this repository is shared more widely.
