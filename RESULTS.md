# Results

## Comparative analysis of mammalian gut microbiome functions associated with species-level cancer mortality

*Status as of 2026-08-02. This document follows the structure `METHODOLOGY.md` Phase 8 specifies, filled in honestly: sections with real findings report them, sections that depend on unrun analysis say so explicitly rather than presenting fabricated numbers. See `PIPELINE_EXECUTION_REPORT.md` for the full phase-by-phase execution log this summarizes.*

---

## Introduction / Hypothesis

Species whose gut microbial community is predicted to produce a higher relative flux of pro-carcinogenic metabolites (secondary bile acids, hydrogen sulfide, colibactin) and/or a lower relative flux of the protective short-chain fatty acid butyrate are hypothesized to have higher species-level cancer mortality risk, after correcting for shared phylogeny and two fixed covariates (litter size, diet category) — see `PREREGISTRATION.md` for the full statement and the rationale for excluding body mass/lifespan as primary covariates.

**This hypothesis has not yet been tested.** No microbiome-feature-vs-cancer-mortality correlation of any kind has been computed. What follows reports what has actually been established (species sample and phylogeny) and is explicit about what hasn't (the statistical test itself).

## Methods (summary — see `METHODOLOGY.md` for full detail)

- **Species overlap (Phase 1):** five published datasets (2 cancer-mortality: Vincze et al. 2022, Boddy et al. 2020; 3 microbiome/metabolome: Milani et al. 2020, Youngblut et al. 2020, Gregor/Guo et al. 2022) were joined on normalized species name, then cross-checked against the GBIF Backbone Taxonomy to catch synonyms a string-based normalizer alone would miss. **Result: 57 species have both cancer-mortality and microbiome data** (up from 55 after the taxonomy pass — see `PIPELINE_EXECUTION_REPORT.md` Phase 1 for the specific synonym pairs found, including one taxonomically contested merge, `Canis rufus`/`Canis lupus`, flagged for manual review rather than accepted uncritically).
- **Phylogeny (Phase 2):** a time-calibrated tree for these 57 species was obtained from TimeTree, after confirming the preferred source (the VertLife/Upham et al. 2019 mammal supertree, which ships as a credible tree set rather than a single point estimate) could not actually be retrieved from this environment (Cloudflare-blocked download, no API token). **56 of 57 species placed on the tree** (`Felis silvestris` excluded — no TimeTree data for it at all); 3 species required a tip relabel to reconcile a genus-level taxonomic revision with this project's data-table naming.
- **Bioinformatics processing, functional annotation, and statistical testing (Phases 3-5):** **not performed.** No raw sequencing data exists in this repository for any of the 56 species, and none of the required tools (Nextflow, fastp, Kraken2, CheckM2, GTDB-Tk, R/phylolm/caper) are installed in this environment. This is a hard data/tooling gap, not a methodological choice.
- **Mechanistic follow-up (Phase 6):** **not performed** — its input (Phase 5 hits) doesn't exist.

## Results

**No Phase 5 statistical results exist.** There is no PGLS effect-size table, no FDR-corrected hit list, and no estimated phylogenetic signal (λ) to report, because the feature matrix Phase 5 requires (per-species microbiome functional data) was never constructed — Phase 3/4, which would have built it, are blocked on missing raw sequencing data.

**No Phase 6 mechanistic predictions exist**, for the same reason.

**What does exist and is reportable:** a validated, phylogenetically-placed, 56-species comparative sample (`data/processed/overlap_species_tree.nwk`, `output/pgls_species_list.csv`) spanning cancer-mortality and microbiome data from five independent published sources — this is the sample size and phylogenetic structure any future PGLS analysis on this project would actually use.

## Limitations

**Sample size / statistical power**, stated here rather than as a discussion afterthought: even once Phase 5 can run, an *n* of 56 species tested against potentially hundreds of functional features has drastically less power than GWAS-style designs achieve with sample sizes in the thousands-to-millions. Realistically only comparatively large effects will be detectable; true-but-modest associations will likely be missed (Type II error) even after Benjamini-Hochberg FDR correction handles the false-positive side.

**Phylogenetic effective sample size:** cannot yet be reported — it depends on the estimated λ from an actual Phase 5 PGLS fit, which hasn't been run. Once it is, λ should be reported alongside the raw species count, since closely related species are not independent data points and 56 raw species overstates the independent information actually available.

**Data provenance caveats:**
- 3 of 5 Phase 1 sources were hand-extracted from manually-downloaded supplementary files (network restrictions in this environment blocked `pmc.ncbi.nlm.nih.gov`, NCBI SRA, and MassIVE).
- The taxonomy-backbone-driven `Canis rufus`→`Canis lupus` merge reflects one side of a genuinely contested taxonomic question and deserves a manual look before being treated as settled.
- The phylogeny is a single TimeTree point estimate, not Upham et al.'s credible tree set, so no sensitivity analysis across phylogenetic uncertainty is currently possible.
- **Phases 3-6 did not run at all** — see `PIPELINE_EXECUTION_REPORT.md` for the specific missing data/tooling. Nothing downstream of Phase 2 in this document reflects a real computed result.

## Next Steps

1. Obtain raw sequencing data (or a validated genome-ID-to-host-species crosswalk for reusing Youngblut et al.'s published MAGs) for the 56 tree species — the single biggest blocker to Phase 3.
2. Get access to a Linux/Docker+Nextflow environment (or cloud compute) with the reference databases CheckM2 and GTDB-Tk need, plus R with `phylolm`/`caper` installed for Phase 5.
3. Once Phase 5 produces real hits, revisit `EXPERIMENTAL_PROPOSAL.md` — its experimental design is already scoped, but its candidate-species section is intentionally left TBD pending real findings.

## Does this support or contradict the Phase 0 hypothesis?

**Neither — it hasn't been tested yet.** Nothing in this pipeline run confirms or refutes the primary hypothesis; the only thing established so far is that a real, taxonomically-vetted, phylogenetically-placed 56-species sample exists to test it on, once Phases 3-6 can actually run.
