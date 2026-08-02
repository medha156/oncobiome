# Preregistration

## Comparative analysis of mammalian gut microbiome functions associated with species-level cancer mortality

**Preregistered on 2026-08-02**, after Phase 1 (data acquisition/species overlap) and Phase 2 (phylogeny) had already run, but **before any Phase 5 statistical correlation test has been attempted.** This ordering is noted explicitly rather than silently glossed over: the ideal preregistration workflow (METHODOLOGY.md Phase 0) commits the hypothesis and covariates before *any* data is touched. Phase 1/2 do not touch the hypothesis test itself — they establish sample size and phylogeny, not cancer-microbiome correlations — so the researcher-degrees-of-freedom risk this document exists to guard against (quietly adjusting the hypothesis or covariates after seeing which one "worked") has not yet been triggered. No microbiome-feature-vs-cancer-mortality correlation of any kind has been computed as of this writing.

---

### 1. Primary hypothesis

> Species whose gut microbial community is predicted (via genome-scale metabolic modeling and flux balance analysis) to produce a higher relative flux of pro-carcinogenic metabolites — secondary bile acids (via microbial 7α-dehydroxylation), hydrogen sulfide, and colibactin — and/or a lower relative flux of the protective short-chain fatty acid butyrate, will have higher species-level cancer mortality risk, after correcting for shared phylogeny and the covariates below.

This is stated at the metabolite-output level, per the project's stated FBA interest (METHODOLOGY.md Phase 0), rather than at the taxonomic-composition or gene-presence level — those are weaker, less mechanistic endpoints available from the same data and are not the primary hypothesis, though they may be reported as secondary/exploratory findings if the metabolite-level analysis (Phase 6) cannot be completed.

### 2. Fixed covariates

Included in every Phase 5 model, decided here rather than after seeing results:

- **Litter size** — from Boddy et al. 2020 [METHODOLOGY.md ref 2]
- **Diet category** (herbivore/carnivore/omnivore) — from Milani et al. 2020 [ref 3]

**Explicitly excluded as a covariate of primary interest:** body mass and maximum lifespan. Both Vincze et al. 2022 and Boddy et al. 2020 tested these against cancer mortality/prevalence and found them non-significant — the empirical basis of Peto's paradox, that large/long-lived animals do not show proportionally higher cancer rates. Since neither source paper found a body-mass/lifespan effect worth correcting for, they are not included as primary covariates here.

**Flagged for Phase 5 sensitivity checks anyway:** body mass and lifespan should still be tested as sensitivity covariates once the real Phase 5 model is run, specifically because this project's ~55-species sample is a different (smaller, non-random) subset of mammals than either source paper's full sample — a null result in a larger sample doesn't guarantee nullity in this specific subsample. This is a check, not a primary-model covariate.

### 3. Planned statistical test

- **Screening pass:** PGLS via R's `phylolm` package across every feature in the Phase 4 feature matrix, response = cancer mortality risk/prevalence, covariates = litter size + diet category, Pagel's λ estimated by maximum likelihood (not fixed). `phylolm` is used for the screening pass specifically for its linear-time scaling with species count (METHODOLOGY.md Phase 5).
- **Multiple-testing correction:** Benjamini-Hochberg false discovery rate across all features tested in the screening pass.
- **Confirmatory pass:** FDR-surviving features only, refit with `caper::pgls()` for full diagnostics (residual plots, profile likelihood for λ).

### 4. Commitment

Preregistered on 2026-08-02, prior to any Phase 5 data pull, feature-matrix construction, or statistical test. The hypothesis, covariate list, and statistical test above will not be changed after Phase 5 results are seen; any post-hoc analysis not specified here will be reported explicitly as exploratory, not confirmatory.
