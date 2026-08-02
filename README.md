# oncobiome

Comparative analysis of mammalian gut microbiome function and species-level
cancer mortality/prevalence (a microbiome-focused look at Peto's paradox:
why large, long-lived animals don't get proportionally more cancer).

## Pipeline overview

| Phase | Script | Output |
|---|---|---|
| 1 | `phase1_data_ingestion.py` | Species-level overlap between cancer-mortality data and microbiome/metabolome data |
| 2 | `phase2_phylogeny.py` | Time-calibrated phylogeny for the overlap species, PGLS-ready |

Both scripts are re-runnable and idempotent — safe to run again as more
source data is added or names are reconciled.

---

## Phase 1 — data ingestion & species overlap

### Data sources (5 papers)

**Cancer-mortality data:**
- Vincze et al. 2022, *Nature* 601:263–267 — cancer risk across 191 mammal species
- Compton/Boddy et al. 2020, *Evolution, Medicine & Public Health* — necropsy-based neoplasia prevalence + life-history traits, 37 species

**Microbiome/metabolome data:**
- Manor et al. 2020, *AEM* 86:e01864-20 — diet/gut-physiology metadata, 77 species
- Youngblut et al. 2020, *mSystems* 5:e01045-20 — metagenome-assembled genomes, 180 species (spans 5 host classes; only ~103 are Mammalia)
- Guo/Gregor et al. 2022, *ISME J* 16:1262–1274 — fecal metabolomics, 25 species

### Choices made

- **3 of 5 sources were hand-extracted, not auto-fetched.** `pmc.ncbi.nlm.nih.gov`, NCBI SRA, and MassIVE were unreachable from this sandbox's network allow-list, so `compton2020_cancer_prevalence.csv`, `youngblut2020_host_species.csv`, and `guo2022_host_species.csv` were extracted from supplementary Excel files downloaded manually and provided directly, rather than fetched programmatically. This is recorded per-source in `DATA_SOURCES[...]["manual_note"]` in the script.
- **Name normalization is string-based only**, not a taxonomy-backbone lookup: underscores → spaces, whitespace collapsed, genus capitalized/epithet lowercased, and trinomial subspecies (e.g. `Equus africanus asinus`) collapsed to their binomial (`Equus africanus`) for joining, since the cancer datasets are species-level only. The full original per-source name is preserved in the output report — only the join key is coarsened.
  - **Why not a taxonomy backbone (NCBI Taxonomy / GBIF) here:** that requires network access and was deferred to Phase 2, where TimeTree's own name resolution ended up serving as a practical substitute (see below) — its "unresolved names" output caught real taxonomic revisions the string normalizer can't.
- **Overlap is defined as the intersection** of species with cancer-mortality data and species with microbiome/metabolome data — this is the number that matters, since it's the actual usable sample size for every downstream statistical step (PGLS, pathway analysis, FBA).
- **The overlap count (55) is explicitly a lower bound**, not a final answer, because of the manual-extraction caveat above and the string-only name matching.
- `data/raw/functional/` (CAZy, KEGG pathway, biosynthetic gene cluster, SGB catalog/phenotype, and annotated metabolite tables from Youngblut/Guo) is **staged but intentionally not wired into Phase 1** — it's the functional-enrichment input for a later analysis phase, once the species sample is locked in.

### Result

- 398 unique species total across all 5 sources
- 209 species have cancer data, 244 have microbiome data
- **55 species have both** (`output/overlap_species_only.csv`, `output/species_overlap_report.csv`)

---

## Phase 2 — time-calibrated phylogeny

### Choice of tree source: TimeTree

Used [TimeTree](http://www.timetree.org) rather than an alternative like Open
Tree of Life or the Upham/VertLife mammal supertree, per direct instruction —
it's the same resource Youngblut et al. 2020 (one of our own source papers)
used, and it returns branch lengths as actual divergence times (millions of
years), which PGLS needs.

### How the tree was actually fetched

TimeTree has **no documented public REST API**. Its "Load a List" → "Prune
Tree" web feature is implemented as two undocumented AJAX endpoints used by
the site's own front-end JavaScript:

```
POST /ajax/prune/load_names/          (multipart file upload, field "file")
POST /ajax/newick/prunetree/download  (form field export=newick)
```

`phase2_phylogeny.py` drives these directly with a `requests` session,
reverse-engineered from TimeTree's public `app.js` — this is exactly what a
browser does when a person uses the "Load a List" tool manually; no
authentication or private access is involved.

**Why this is flagged as a real caveat, not just an implementation detail:**
these are unversioned, undocumented endpoints TimeTree could change at any
time without notice, which would silently break this script. If it stops
working, the fallback is to use the interactive "Load a List" tool by hand
and drop the resulting Newick file at the path `RAW_TREE_PATH` points to.

### Why the tree cross-check was run before any statistics

Submitting the species list to TimeTree doubles as a taxonomic sanity check —
TimeTree's own name resolution surfaces cross-dataset naming problems the
Phase 1 string normalizer can't catch (synonyms, genus reassignments,
taxa TimeTree has no data for at all). Running this early, before PGLS,
means these get caught and decided on deliberately instead of silently
corrupting a downstream join.

### What the cross-check found, and the choices made in response

Of the 55 Phase 1 overlap species, 54 were placed on the tree:

| Species | Issue | Choice made | Why |
|---|---|---|---|
| `Felis silvestris` | TimeTree: "insufficient data ... to place this taxon" — dropped entirely, no tip at all | **Excluded** from the final PGLS species list | No tip exists to graft or rename. The alternative — manually grafting a tip at a literature-sourced divergence time — would inject an assumption not backed by TimeTree's own data, so it wasn't done silently or by default. |
| `Callithrix pygmaea` | Tree tip is `Cebuella_pygmaea` (real genus reassignment, pygmy marmoset), **not flagged** by TimeTree's own unresolved-names list | **Tip relabeled** to `Callithrix_pygmaea` in the final tree | PGLS tools (e.g. R's `ape`/`caper`) join tree tips to trait-table rows by exact string match. The Phase 1 data tables already use `Callithrix pygmaea` as the join key, so the tip was renamed to match the existing pipeline's naming rather than renaming the naming convention everywhere else. Topology/branch length is untouched — this is a join-key fix only, logged explicitly in `output/pgls_species_list.csv`'s `note` column so it stays auditable, not hidden. |
| `Macropus eugenii` | Tree tip is `Notamacropus_eugenii` (genus split), also **not flagged** by TimeTree | **Tip relabeled** to `Macropus_eugenii` | Same reasoning as above. |
| `Gazella subgutturosa` | TimeTree flagged this itself: substituted with `Gazella dorcas`'s branch as a data proxy | **Kept, not excluded** — but the substitution is carried into the `note` column | TimeTree explicitly told us about this one (unlike the two silent renames above), so no detective work was needed — just surfacing the caveat downstream: this species' branch length is a congener's data, not species-specific, which matters for how much weight to put on it in PGLS. |
| `Giraffa camelopardalis` | Same situation: substituted with `Giraffa reticulata`'s branch | **Kept, with note carried forward** | Same reasoning. |

**Net result: 54/55 species usable for PGLS**, with 2 label fixes applied and
1 species excluded — all decisions logged, not silently absorbed.

### Outputs

- `data/raw/phylogeny/timetree_prunetree_raw.nwk` — raw tree exactly as returned by TimeTree, untouched
- `data/raw/phylogeny/timetree_pgls_ready.nwk` — same tree, tips relabeled per the table above; this is the one to use for PGLS
- `data/raw/phylogeny/timetree_submitted_species_list.txt` — snapshot of exactly what was submitted, for reproducibility
- `output/phylogeny_species_match_report.csv` — full per-species tree-placement detail
- `output/pgls_species_list.csv` — final 54-species list with sources and caveat notes, ready to join against the Phase 1 trait tables
- `logs/phase2_*.log` — full timestamped run log

---

## Known limitations / next steps

- Phase 1's overlap is a lower bound (manual-extraction sources, string-only name matching); wiring a proper taxonomy backbone (NCBI Taxonomy/GBIF) instead of relying on TimeTree's incidental name resolution would tighten this.
- The TimeTree fetch depends on undocumented endpoints and could break if TimeTree changes their site (see Phase 2 caveat above).
- `data/raw/functional/` (CAZy, KEGG, BGC, SGB, metabolite tables) is staged but not yet used — planned for a functional-enrichment phase once the 54-species PGLS sample above is finalized.
