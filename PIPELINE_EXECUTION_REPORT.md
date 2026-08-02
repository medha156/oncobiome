# Pipeline Execution Report

## What happened when `PIPELINE_PROMPTS.md` was run end-to-end, Phase 0 through Phase 8

**Run date: 2026-08-02.** This document is the log of every decision and change made while executing the prompts in `PIPELINE_PROMPTS.md` in order. Per that file's own instructions, phases that turned out to be genuinely blocked (missing data or missing tooling) are reported as blocked rather than faked — nothing in this report or in the repo's output files reflects a fabricated result.

**Bottom line up front:** Phases 0, 1, and 2 executed for real and produced real, verifiable outputs. Phases 3 through 6 are **blocked** — not on a fixable bug, but on two hard constraints of this environment: there is no raw sequencing data (FASTQ/SRA) anywhere in this repo for any species, and none of the required bioinformatics tools (Nextflow, fastp, Kraken2, CheckM2, GTDB-Tk, CarveMe, gapseq, COBRApy, MICOM, R) are installed. Phases 7 and 8 were written up as far as they honestly could be without fabricating scientific findings that don't exist yet.

---

## Phase 0 — Scope the question

**Status: Done.** `PREREGISTRATION.md` written.

**Decisions made:**
- Hypothesis fixed at the metabolite-output level (secondary bile acids, H2S, colibactin vs. butyrate), per the project's stated FBA interest, not the weaker taxonomic-composition level.
- Covariates fixed: litter size (Boddy et al. 2020), diet category (Milani et al. 2020). Body mass/lifespan explicitly excluded as primary covariates (both source papers found them non-significant — Peto's paradox), but flagged for mandatory sensitivity-check re-testing in Phase 5, since this project's ~55-species subsample isn't guaranteed to reproduce a null result found in a larger sample.
- **Chronology caveat, stated explicitly in the document itself:** this preregistration was written *after* Phase 1 and Phase 2 had already executed in earlier work on this repo, which technically breaks the "commit before touching any data" ideal `METHODOLOGY.md` Phase 0 describes. This is flagged rather than hidden: Phase 1/2 don't touch the actual hypothesis test (no cancer-vs-microbiome correlation has been computed at any point), so the degrees-of-freedom risk preregistration guards against hasn't actually been exploited — but the ordering itself is a real deviation from ideal practice, worth knowing about.

---

## Phase 1 — Data acquisition and species overlap check (re-run + taxonomic backbone resolution)

**Status: Done.** Step 1 (naming fix) was already completed in prior work. Step 2 (taxonomic backbone resolution) executed this run via a new script, `phase1b_taxonomy_resolution.py`.

**Feasibility check performed:** both GBIF Backbone Taxonomy (`api.gbif.org`) and NCBI Taxonomy (`eutils.ncbi.nlm.nih.gov`) were reachable from this environment (confirmed via direct HTTP request, HTTP 200 on both). **GBIF was used** — its `/species/match` endpoint returns the resolved accepted binomial in a single call per species, versus NCBI E-utilities' two-step search+fetch flow, which is simpler to drive at the ~340-species scale this step needed.

**What it did:** every species Phase 1 found in only *one* of the two data categories (154 cancer-only, 189 microbiome-only — 343 total) was resolved against GBIF. Where a cancer-only species and a microbiome-only species resolved to the *same* accepted binomial under different spellings, that's a real overlap species the string-based normalizer missed.

**Result — before/after counts:**
| Stage | Overlap count |
|---|---|
| Original Phase 1 run (string matching only) | 55 |
| After the Milani/Boddy naming fix (unchanged, as expected) | 55 |
| **After GBIF taxonomic backbone resolution** | **57** |

**2 new overlap species found:**
- **`Canis lupus`** — cancer-mortality data had it as `Canis rufus` (red wolf); GBIF backbone resolves `Canis rufus` to `Canis lupus` (SYNONYM, 98% confidence, EXACT match). **Flagged as taxonomically contested, not a clean win**: the red wolf's status as a distinct species vs. a wolf/coyote hybrid lineage is an active, genuinely disputed question in the literature (IUCN/USFWS treat it as distinct in conservation contexts; some genomic studies argue for a hybrid origin folding it into *C. lupus*). This merge should get a manual second look before being treated as settled, rather than accepted purely because one taxonomy backbone says so.
- **`Taurotragus oryx`** — cancer data had `Tragelaphus oryx` (common eland); GBIF resolves this to `Taurotragus oryx` (SYNONYM, 98% confidence, EXACT match) — a standard, uncontroversial genus reassignment, same pattern as the `Cebuella`/`Callithrix` and `Notamacropus`/`Macropus` renames Phase 2 already caught independently.

**4 additional cases logged** (don't change the overlap count, but widen source support for species already counted): `Cervus canadensis`→`Cervus elaphus`, `Dama mesopotamica`→`Dama dama`, `Equus burchelli`→`Equus quagga`, `Lama guanicoe`→`Lama glama`.

**17 species had no usable GBIF match at all** — left as unmatched, not silently dropped (see `output/phase1b_taxonomy_resolution_report.csv` for the full per-species audit trail, including every name that *didn't* resolve).

**Outputs updated:** `output/species_overlap_report.csv`, `output/overlap_species_only.csv` (now 57 species), plus the new `output/phase1b_taxonomy_resolution_report.csv` audit trail and `logs/phase1b_*.log`.

---

## Phase 2 — Build a species-level phylogenetic tree

**Status: Done**, re-run against the updated 57-species list.

**VertLife/Upham feasibility check (done before defaulting to TimeTree, per the prompt's explicit instruction):**
- The full Dryad data package for Upham et al. 2019 (`doi:10.5061/dryad.tb03d03`) is **5.5 GB** — impractical to bulk-download just to extract one 57-species subset.
- The smaller, specifically relevant file (`Data_S4_patchClade_results_and_MCC.zip`, ~4.5 MB, the maximum-clade-credibility tree) **could not be downloaded**: Dryad's REST API file-download endpoint requires an OAuth bearer token this environment doesn't have (401 Unauthorized), and the public web download link is behind a Cloudflare JavaScript challenge ("Validating..." page) that a plain HTTP client cannot pass.
- The interactive VertLife species-subsetting tool (`vertlife.org/phylosubsets`) does load as a page, but isn't a simple stateless request/response flow the way TimeTree's "Load a List" AJAX endpoints turned out to be — no quick scriptable path was found there either.
- **Decision: fell back to TimeTree**, as `METHODOLOGY.md` explicitly instructs when VertLife isn't reachable. **This is logged as a real limitation, not a stylistic choice**: it means Phase 5 (whenever it can run) will use a single point-estimate tree, not Upham et al.'s credible tree set, and cannot re-run PGLS across tree replicates to quantify sensitivity to phylogenetic uncertainty.

**Cross-check results on the 57-species list:** 56 placed on the tree, same pattern as before plus one more:
- `Felis silvestris` — still no TimeTree data at all, excluded.
- `Callithrix pygmaea` → tree tip `Cebuella_pygmaea` (silent genus rename, tip relabeled to match).
- `Macropus eugenii` → tree tip `Notamacropus_eugenii` (silent genus rename, tip relabeled to match).
- **New: `Taurotragus oryx`** → tree tip `Tragelaphus_oryx`. Interestingly, **TimeTree and GBIF disagree** on which genus name is currently accepted for this species (TimeTree's tip uses the old `Tragelaphus oryx`, GBIF's backbone says `Taurotragus oryx` is accepted) — a real example of two legitimate taxonomic authorities not agreeing, not a data error. Tip relabeled to `Taurotragus_oryx` to match this project's Phase 1 join key.
- `Gazella subgutturosa` / `Giraffa camelopardalis` — same TimeTree-flagged proxy substitutions as before (branch length borrowed from `Gazella dorcas` / `Giraffa reticulata`).

**Net: 56/57 species usable for PGLS.**

**Outputs:** `data/raw/phylogeny/timetree_pgls_ready.nwk` (as before) **plus a new canonical path, `data/processed/overlap_species_tree.nwk`**, matching the path later phases (Phase 5's prompt) expect. `output/pgls_species_list.csv` updated to 56 species.

---

## Phase 3 — Bioinformatics processing of the microbiome data

**Status: BLOCKED. Feasibility checked first, nothing was run that would fail partway through.**

Checked in this environment:
| Requirement | Status |
|---|---|
| Raw sequencing data (FASTQ/SRA) for any of the 56 tree species | **None found anywhere in the repo** (`data/raw/` contains only species-level summary CSVs — taxonomy, counts, gene-family enrichment statistics — never raw reads) |
| Nextflow | Not installed |
| Docker | Client installed (v29.2.1), but **daemon not running** (`failed to connect to the docker API ... daemon is running?`) |
| fastp, Kraken2, CheckM2, GTDB-Tk | None installed |
| Disk space | 93 GB free — not the constraint |

**Reuse-existing-MAGs path also checked and also blocked:** Youngblut et al. 2020's published genome IDs (`data/raw/functional/youngblut_sgb_catalog.csv`) don't consistently encode host species. Some genome-ID prefixes do embed a recognizable species name (e.g. `Cebus_capucinus_imitator_PRJNA485217_SSR047`), but most don't (e.g. `Ash_p1`, `Bissell_p2`, `Cam_p1`) — those look like individual-animal codes that would need the paper's original sample-metadata crosswalk (linking genome/MAG ID → sample ID → host species) to resolve, and that crosswalk isn't among the files currently staged in this repo. Building that mapping by guessing from filename patterns alone would risk silently mis-assigning genomes to the wrong host species, so it wasn't attempted.

**Conclusion: Phase 3 cannot run on real data in this environment**, regardless of which tool is used, because there is nothing to feed it. This is a data-availability gap, not a tool-installation gap that pip/conda could fix.

**What would unblock this:** either (a) raw FASTQ/SRA accessions for the tree species, or (b) the Youngblut et al. supplementary genome-ID-to-host-species crosswalk table, plus a Linux/Docker+Nextflow environment (or a cloud compute environment) to actually run `nf-core/mag`, CheckM2, and GTDB-Tk.

---

## Phase 4 — Functional annotation and feature matrix construction

**Status: BLOCKED**, for the same root cause as Phase 3: no per-species MAGs exist to annotate (Phase 4's entire input is Phase 3's output).

**Tool availability checked anyway, for completeness:** no eggNOG-mapper, no KofamKOALA, no dbCAN3 in this environment.

**Considered and rejected as a substitute:** the already-staged Youngblut enrichment tables (`data/raw/functional/youngblut_*_enriched_*.csv`, KEGG/CAZy tables) are *group-level* enrichment test results (host-vs-environment, mammal-vs-other), not a per-species feature table — they answer "is this CAZy family enriched in hosts vs. environment," not "how much of this CAZy family does species X's community carry." Using them as a stand-in for the Phase 4 deliverable would misrepresent a two-group comparison as species-level data, so this wasn't done.

**Conclusion:** no `data/processed/species_feature_matrix.csv` exists, and none was fabricated.

---

## Phase 5 — Statistical correlation with phylogenetic correction

**Status: BLOCKED** — its entire input (`data/processed/species_feature_matrix.csv`) doesn't exist.

**Tool availability checked anyway:** no R, no Rscript in this environment (needed for `phylolm`/`caper`). Even if a real feature matrix existed, the statistical tool stack for the exact method `METHODOLOGY.md` specifies isn't currently installed here either.

**No PGLS was run. No p-values, effect sizes, or FDR results exist**, and none were invented. `output/phase5_pgls_results.csv` and `output/phase5_significant_hits.csv` were not created.

---

## Phase 6 — Mechanistic follow-up on the hits

**Status: BLOCKED** — its input (`output/phase5_significant_hits.csv`) doesn't exist, because Phase 5 didn't run.

**Tool availability checked anyway:** `cobra` (COBRApy), `carveme`, and `micom` Python packages are not installed; no `gapseq` binary either.

**No genome-scale metabolic models, no FBA predictions, and no "species X is predicted to overproduce compound Y" claims exist.** `output/phase6_fba_predictions.csv` was not created, and no mechanistic claims were fabricated.

---

## Phase 7 — Experimental validation proposal

**Status: Partially written.** `EXPERIMENTAL_PROPOSAL.md` was created, but its candidate-species/metabolite section is explicitly marked **TBD** rather than filled in with invented names. The general experimental *design* (germ-free FMT vs. alternatives, sample-size/pre-registration discipline responding to Walter et al. 2020, ARRIVE 2.0 commitment) doesn't depend on which species turn out to be hits, so that part was written for real. The part that names 1-3 specific candidate species can only be filled in once Phase 5/6 produce real hits — doing so now would mean inventing a scientific finding, which this report is explicit about not doing.

---

## Phase 8 — Write-up

**Status: Partially written.** `RESULTS.md` was created following the Phase 8 structure (Introduction/Hypothesis, Methods, Results, Limitations, Next Steps), but its Results section honestly states that no Phase 5/6 findings exist yet, rather than presenting fabricated statistics. The Methods section accurately summarizes what Phases 1-2 actually established (57-species overlap, 56-species pruned tree). The Limitations section covers both the standard sample-size/power framing `METHODOLOGY.md` calls for *and* the specific blocked-phases situation documented in this report.

---

## Overall consistency check

- `output/overlap_species_only.csv`: 57 species (verified: 58 lines including header).
- `output/pgls_species_list.csv`: 56 species (Felis silvestris excluded, matches Phase 2's tree placement exactly).
- `data/processed/overlap_species_tree.nwk` and `data/raw/phylogeny/timetree_pgls_ready.nwk`: identical content, both reflect the 56-species tree with the same 3 tip relabels (Callithrix pygmaea, Macropus eugenii, Taurotragus oryx).
- No file in this repo claims a PGLS result, an FBA prediction, or a specific experimental candidate that wasn't actually produced by a real, logged computation.

## What's needed to unblock Phases 3-6

1. Raw sequencing data (FASTQ/SRA accessions) for the 56 tree species not already covered by clean, species-mapped published MAGs — this is the single biggest gap.
2. A Linux environment with Docker/Nextflow actually running (or a cloud compute environment), plus the reference databases CheckM2/GTDB-Tk need.
3. R with `phylolm` and `caper` installed, for whenever a real feature matrix exists.
4. `cobra`, `carveme`, `micom` (pip-installable) and `gapseq` for Phase 6, once Phase 5 produces real hits to follow up on.
