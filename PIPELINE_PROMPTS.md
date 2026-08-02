# Execution Prompts

## A ready-to-paste prompt for each phase of the oncobiome pipeline

These prompts turn `METHODOLOGY.md` into something directly executable. Each one is meant to be handed to Claude (in a fresh conversation, or in sequence in this same one) at the point where the prior phase's deliverable already exists in the repo. They assume the `oncobiome` repository structure already in place: `phase1_data_ingestion.py`, `data/raw/`, `data/raw/functional/`, and `output/species_overlap_report.csv` / `output/overlap_species_only.csv`.

A few notes before using these:

- **Run them in order.** Each prompt explicitly names the prior phase's output file(s) as its input, so skipping ahead will just produce an error or a prompt for missing data.
- **Some phases need compute this sandbox doesn't have.** Phase 3 (assembly/binning) and Phase 6 (GEM reconstruction/FBA) call for tools like `nf-core/mag`, GTDB-Tk, and CarveMe that need either a large reference database download, a Nextflow/Docker environment, or both — more than a lightweight sandbox typically provides. Each of those prompts includes a fallback instruction to first check feasibility and report back rather than silently failing or fabricating output.
- **Naming issue: FIXED 2026-08-02.** `phase1_data_ingestion.py` previously had `manor2020_microbiome` and `compton2020_cancer` as misnomers (real authors: Milani et al. and Boddy et al., respectively — see `METHODOLOGY.md`'s closing note). These have been renamed to `milani2020_microbiome` and `boddy2020_cancer` (dict keys, citations, and underlying file names `milani2020_host_species.csv` / `boddy2020_cancer_prevalence.csv`), and Phase 1 + Phase 2 were re-run to confirm the overlap count is unchanged (still 55). The Phase 1 prompt below's step 1 (the rename) is therefore already done — only step 2 (taxonomic backbone resolution pass) remains outstanding.
- **Adjust species/sample counts as the real pipeline runs.** The counts referenced (55-species overlap, 103 mammal species from Youngblut et al., etc.) reflect the state of the repo as of this writing; update them if Phase 1 is re-run and the overlap changes.

---

## Phase 0 — Scope the question

```
Read METHODOLOGY.md in this repository, specifically the Phase 0 section. Based on
it, draft a short preregistration document for this project and save it as
PREREGISTRATION.md in the repo root.

Include:
1. The primary hypothesis, stated as a single testable sentence, at the
   metabolite-output level (per the project's FBA interest) — e.g., framed around
   specific microbial metabolite classes (secondary bile acids, H2S, colibactin,
   SCFAs) and species-level cancer mortality risk.
2. The fixed covariate list: litter size (from Boddy et al. 2020) and diet category
   (from Milani et al. 2020), and note explicitly that body mass/lifespan will NOT be
   included as a covariate of interest since Vincze et al. 2022 and Boddy et al. 2020
   both found these non-significant (consistent with Peto's paradox) — but flag that
   they should still be checked as sensitivity covariates in Phase 5.
3. The planned statistical test (PGLS via phylolm for screening, caper for
   confirmatory models — see METHODOLOGY.md Phase 5) and the multiple-testing
   correction to be used (Benjamini-Hochberg FDR).
4. A dated commitment line ("Preregistered on [date], prior to Phase 1 data
   pulls/re-runs") — use today's date.

Keep it to about one page. This is a lightweight, repo-committed preregistration,
not a formal OSF submission (see METHODOLOGY.md Phase 0 "Alternatives considered"
for why). When done, tell me the hypothesis statement you settled on so I can
confirm before we move to Phase 1.
```

---

## Phase 1 — Data acquisition and species overlap check (re-run + cleanup)

```
This repo already has a working Phase 1 pipeline (phase1_data_ingestion.py) that
found a 55-species overlap across five source datasets. Step 1 below (the naming
fix) is already done as of 2026-08-02 — DATA_SOURCES keys are now
`milani2020_microbiome` and `boddy2020_cancer`, matching their real authors
(Milani et al. 2020 AEM, Boddy et al. 2020 EMPH), and the underlying files are
`milani2020_host_species.csv` / `boddy2020_cancer_prevalence.csv`. The overlap
count was confirmed unchanged (55) after the rename. Only step 2 remains:

2. Add a second-pass taxonomic name resolution step, per METHODOLOGY.md Phase 1
   "Alternatives considered": after the existing string-based normalizer runs, take
   any species names that failed to find a match across sources and attempt to
   resolve them via a taxonomic backbone (GBIF Backbone Taxonomy API or NCBI
   Taxonomy, whichever is reachable from this environment — check network access
   first and report which one works). Log every name that gets resolved this way
   (before name -> after name -> source dataset), and log any names that still fail
   to resolve. Add these newly-resolved matches to the overlap if they weren't
   already counted.

Update output/species_overlap_report.csv and output/overlap_species_only.csv with
the results. Give me detailed before/after counts: original overlap (55), overlap
after the naming fix (should be unchanged), and overlap after the taxonomic
backbone pass (may increase). Flag anything that looks wrong rather than assuming
success.
```

---

## Phase 2 — Build a species-level phylogenetic tree

```
Using output/overlap_species_only.csv from Phase 1 as the target species list,
build a pruned phylogenetic tree for this project, following METHODOLOGY.md
Phase 2.

1. Obtain a mammal phylogeny covering the overlap species. Preferred source: the
   VertLife/Upham et al. 2019 mammal supertree (credible tree set). Check whether
   this is reachable/downloadable from this environment; if not, fall back to
   TimeTree.org and note explicitly in your output that you used the fallback and
   why (per METHODOLOGY.md's "Alternatives considered" for Phase 2, this loses the
   multi-tree uncertainty quantification Upham et al. provides — flag this as a
   limitation rather than silently substituting).

2. Prune the tree to exactly the overlap species list using something equivalent to
   R's `ape::drop.tip()` (or a Python equivalent if R isn't available in this
   environment — check first, tell me which you're using).

3. In both directions: drop any tree tips not in the overlap list, AND flag any
   overlap species not found in the tree. For flagged species, try matching against
   tree tip labels using the same taxonomic-backbone resolution approach from the
   Phase 1 prompt (this is the second independent QC pass METHODOLOGY.md Phase 2
   describes) before giving up and logging them as unmatched.

4. Save the final pruned tree as data/processed/overlap_species_tree.nwk (Newick
   format), and write a short log (logs/phase2_<timestamp>.log) recording: starting
   species count, tips dropped, species unmatched to any tree tip (with reasons),
   and final tree species count.

Report back: final tree species count vs. the Phase 1 overlap count, and list any
species that didn't make it into the tree.
```

---

## Phase 3 — Bioinformatics processing of the microbiome data

```
Before doing any heavy compute, check feasibility first: this phase normally uses
nf-core/mag (Nextflow), CheckM2, and GTDB-Tk, each of which needs either large
reference databases (tens of GB) or a container runtime. Check what's actually
available in this environment (Nextflow, Docker/Singularity, disk space, network
access to database downloads) and report back BEFORE attempting to run anything
that will fail partway through. If the full pipeline isn't runnable here, say so
clearly and propose what can be done instead (e.g., processing a single small
representative sample, or preparing the pipeline config/scripts for a person to
run on adequate compute elsewhere).

Assuming raw sequencing data or reused MAGs are available for at least a subset of
the Phase 2 tree's species, run the following per METHODOLOGY.md Phase 3:

1. First, check whether any of the target species are already covered by
   Youngblut et al. 2020's published MAGs (see data/raw/functional/
   youngblut_sgb_catalog.csv, 1,522 SGBs from 180 host species). For any overlap
   species covered there, treat this step as done — do not re-assemble — and note
   which species these are.

2. For species requiring new processing: run fastp for adapter/quality trimming,
   then host-read removal against each sample's own host reference genome (or
   nearest available relative — look this up per species and log which reference
   was used for which sample).

3. Run Kraken2 as a fast taxonomic sanity check on QC'd reads (not MetaPhlAn4 — see
   METHODOLOGY.md Phase 3 "Alternatives considered" for why Kraken2 is preferred at
   this stage). Flag any sample whose top hits don't match the expected host
   species' known gut flora as a possible swap/contamination issue.

4. Assemble and bin via nf-core/mag if feasible (per the feasibility check above).
   Quality-filter resulting bins with CheckM2 (completeness/contamination), and
   assign taxonomy with GTDB-Tk.

Log everything — per-sample QC stats, taxonomic sanity-check results, CheckM2
scores, GTDB-Tk assignments — the same way phase1_data_ingestion.py logs its
steps (dual file+console handlers, INFO console / DEBUG file). Save quality-passed
MAGs to data/processed/mags/. Give me an honest status report at the end: what
ran, what didn't, and what's blocked on compute/data availability.
```

---

## Phase 4 — Functional annotation and feature matrix construction

```
Using the quality-filtered MAGs from Phase 3 (data/processed/mags/, plus reused
Youngblut et al. MAGs/SGBs where applicable — see
data/raw/functional/youngblut_sgb_phenotypes.csv and the KEGG/CAZy enrichment
tables already in data/raw/functional/), build the Phase 4 feature matrix per
METHODOLOGY.md.

1. Annotate each MAG with both eggNOG-mapper and KofamKOALA (check tool
   availability in this environment first and report). Cross-validate hits
   specifically for the biologically important gene set: secondary bile acid
   7-alpha-dehydroxylase genes (bai operon), butyrate synthesis genes (but/buk),
   and colibactin-associated pks genes. Flag any disagreement between the two
   annotation tools for these specific genes.

2. Run dbCAN3 separately for CAZyme family and substrate-level prediction.

3. Collapse gene-level annotations to a per-species feature matrix: pathway
   presence/absence, gene family relative abundance, and specific counts for the
   bile-acid/butyrate/colibactin genes above, one row per Phase 2 tree species.

4. Apply centered log-ratio (CLR) transformation to the resulting abundance
   features — NOT rarefaction (see METHODOLOGY.md Phase 4 "Alternatives
   considered" for why). Use a standard pseudocount approach for zero-handling
   before the log transform, and document what pseudocount you used.

Save the final table as data/processed/species_feature_matrix.csv (species as
rows, CLR-transformed features as columns, plus a few raw/untransformed columns
for the bile-acid/butyrate/colibactin gene counts specifically, kept both raw and
transformed since Phase 6 will want raw counts too). Log row/column counts, how
many species had at least one MAG available, and any species that ended up with
no functional data at all (and therefore need to be dropped from Phase 5).
```

---

## Phase 5 — Statistical correlation with phylogenetic correction

```
Using data/processed/species_feature_matrix.csv (Phase 4), the pruned tree
data/processed/overlap_species_tree.nwk (Phase 2), the cancer mortality data from
Phase 1's sources, and the covariates fixed in PREREGISTRATION.md (Phase 0), run
the Phase 5 statistical analysis per METHODOLOGY.md.

1. Screening pass: fit PGLS via phylolm (R) across every feature column in the
   feature matrix, with cancer mortality risk (or prevalence, depending on which
   source covers a given species — log which metric was used per species) as the
   response, litter size and diet category as covariates, and lambda estimated by
   maximum likelihood rather than fixed. Use phylolm specifically for this pass
   because of its linear-time scaling (see METHODOLOGY.md Phase 5 "Alternatives
   considered").

2. Apply Benjamini-Hochberg FDR correction across all features tested in the
   screening pass. Report the number of features surviving at FDR < 0.05 (or
   whatever threshold you and I agree makes sense given the ~55-species n — flag
   if this feels underpowered and say so explicitly rather than just reporting a
   number).

3. Confirmatory pass: for FDR-surviving features only, refit using caper's pgls()
   for full diagnostics (residual plots, profile likelihood for lambda). Save
   diagnostic plots to output/phase5_diagnostics/.

4. Flag phylogenetic outliers: species whose residuals are notably large relative
   to the fitted model, for both the full feature set and specifically for the
   FDR-surviving hits. List these with a brief note on each (known unusual
   biology, e.g. naked mole-rat-style cancer resistance, vs. possible data quality
   issue).

Save the full results table (feature, effect size, t-value, raw p-value,
FDR-adjusted p-value, estimated lambda) to output/phase5_pgls_results.csv, and a
separate short list of just the FDR-surviving hits to
output/phase5_significant_hits.csv. Report the hit list, effect sizes, and
outlier species directly in your response, not just in the saved files.
```

---

## Phase 6 — Mechanistic follow-up on the hits

```
For the species/features in output/phase5_significant_hits.csv (Phase 5), run the
Phase 6 mechanistic follow-up per METHODOLOGY.md. Check feasibility first (CarveMe,
gapseq, COBRApy, and MICOM all need to be installed, and gapseq/MICOM in particular
can be slow — report what's available before running anything at scale).

1. Build genome-scale metabolic models (GEMs) with CarveMe for the MAGs belonging
   to the Phase-5 hit species (bulk reconstruction, including CarveMe's
   community-model mode). For the specific bile-acid/butyrate/colibactin pathway
   regions identified in Phase 4, do a second, targeted rebuild with gapseq for
   higher per-pathway accuracy (see METHODOLOGY.md Phase 6 "Alternatives
   considered" for why both tools are used, each for a different purpose).

2. Run community-level FBA with MICOM (not single-organism FBA, and not
   SteadyCom — see METHODOLOGY.md Phase 6 for why), constraining the nutrient
   input to each species' actual diet category from the Phase 1 microbiome
   datasets (herbivore/carnivore/omnivore), not a generic default medium.

3. Cross-reference the predicted flux output against the curated cancer-relevant
   metabolite list from METHODOLOGY.md Phase 4/6: secondary bile acids (via
   7-alpha-dehydroxylation), hydrogen sulfide, colibactin, and butyrate/other
   SCFAs. For each Phase-5 hit species, report whether the FBA output shows
   elevated predicted flux through any of these specific compounds.

Save GEMs to data/processed/gems/, FBA results to
output/phase6_fba_predictions.csv (species, diet constraint used, predicted flux
per compound of interest), and give me a short list of specific,
falsifiable mechanistic claims this phase generated — e.g. "species X's gut
community is predicted to overproduce deoxycholic acid, consistent with its
Phase 5 cancer-mortality association." These are the candidates Phase 7 would
validate experimentally.
```

---

## Phase 7 — Experimental validation proposal (write-up only, not wet-lab work)

```
This phase isn't something to execute computationally — it's a proposal document,
per METHODOLOGY.md Phase 7. Using output/phase6_fba_predictions.csv and the
mechanistic claims from Phase 6, draft an experimental validation proposal and
save it as EXPERIMENTAL_PROPOSAL.md in the repo root.

Include:
1. Which 1-3 candidate species/metabolite claims from Phase 6 are being proposed
   for validation, and why these were prioritized over the rest of the hit list
   (clearest mechanistic story per METHODOLOGY.md Phase 6/7 — state the specific
   reasoning for each).
2. The proposed design: FMT into germ-free (gnotobiotic) mice, comparing
   transplants from high- vs. low-cancer-mortality-associated source microbiomes,
   per METHODOLOGY.md Phase 7. Explicitly state why germ-free FMT was chosen over
   in vitro bioreactor systems and antibiotic-depletion mouse models (see
   METHODOLOGY.md's "Alternatives considered" for this phase), and mention the
   bioreactor pre-screen as an optional cheaper first step.
3. A pre-specified phenotype/endpoint and a sample-size justification, written
   with Walter et al. 2020's critique of this exact experimental design in mind
   (95% of published HMA-rodent studies report phenotype transfer, a rate that's
   probably inflated by publication bias and underpowering — the proposal should
   explicitly avoid repeating that pattern).
4. A commitment to ARRIVE 2.0 reporting from the design stage, not retrofitted
   later.

Keep this to 1-2 pages. This is a scoping document for a follow-on study, not
something we're running in this pipeline — make that framing explicit in the
document itself.
```

---

## Phase 8 — Write-up

```
Using everything produced so far — PREREGISTRATION.md (Phase 0),
output/species_overlap_report.csv (Phase 1), data/processed/overlap_species_tree.nwk
(Phase 2), Phase 3's processing logs, data/processed/species_feature_matrix.csv
(Phase 4), output/phase5_pgls_results.csv and output/phase5_significant_hits.csv
(Phase 5), output/phase6_fba_predictions.csv (Phase 6), and
EXPERIMENTAL_PROPOSAL.md (Phase 7) — draft the project's results write-up per
METHODOLOGY.md Phase 8, and save it as RESULTS.md in the repo root.

Follow METHODOLOGY.md Phase 8's reporting requirements exactly:
1. Report effect sizes and confidence intervals for every Phase 5 tested feature,
   not just p-values or FDR status.
2. State the sample-size/power limitation explicitly, in the methods/introduction,
   not as a discussion afterthought — be specific about the species count, the
   number of features tested, and how this compares in power to GWAS-style designs
   (see METHODOLOGY.md Phase 8 for the specific framing).
3. Report the phylogenetic effective sample size (estimated lambda from Phase 5)
   alongside the raw species count.
4. Keep Phase 5 statistical associations, Phase 6 mechanistic/FBA predictions, and
   the Phase 7 proposal clearly and separately labeled — do not let causal language
   drift upward between sections. A PGLS hit is not the same claim as an FBA
   prediction, which is not the same claim as validated causality.

Structure it as: Introduction/Hypothesis (from Phase 0), Methods (summarizing
Phases 1-4, linking out to METHODOLOGY.md for full detail rather than repeating
it), Results (Phase 5 + Phase 6), Limitations (sample size/power, phylogenetic
effective n, data provenance caveats), and Next Steps (pointing to
EXPERIMENTAL_PROPOSAL.md). When done, tell me if any part of the results actually
supports or contradicts the original Phase 0 hypothesis, plainly, in a sentence or
two.
```
