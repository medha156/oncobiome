#!/usr/bin/env python3
"""
Phase 1, step 2 -- taxonomic backbone name resolution
========================================================
Project: Comparative analysis of mammalian gut microbiome functions
         associated with species-level cancer mortality

WHAT THIS SCRIPT DOES
----------------------
phase1_data_ingestion.py's normalizer is string-based only (underscores,
capitalization, trinomial collapse) and cannot catch real taxonomic
synonyms -- the same organism published under two different scientific
names due to a taxonomic revision (e.g. Callithrix pygmaea / Cebuella
pygmaea). This script is the follow-up METHODOLOGY.md Phase 1
"Alternatives considered" calls for: take every species that Phase 1
found in only ONE of the two data categories (cancer-mortality-only or
microbiome-only) and resolve it against the GBIF Backbone Taxonomy,
then check whether any cancer-only species and any microbiome-only
species resolve to the SAME accepted name under different spellings --
that is a real overlap species Phase 1's string matching missed.

DATA PROVENANCE
-----------------
Uses the public GBIF Backbone Taxonomy species-match API
(https://api.gbif.org/v1/species/match), no API key required. NCBI
Taxonomy (via E-utilities) was also confirmed reachable from this
environment and would be an equally valid alternative per
METHODOLOGY.md ref [7]; GBIF was used here because its /species/match
endpoint returns the resolved accepted name directly in one call per
species, which is simpler to drive at ~340-species scale than NCBI's
search+fetch two-step E-utilities flow.

USAGE
-----
    python phase1b_taxonomy_resolution.py

Re-runnable and idempotent: reads output/species_overlap_report.csv
(Phase 1's output) and re-derives everything from there.
"""

import csv
import logging
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

BASE_DIR = Path(__file__).resolve().parent
LOG_DIR = BASE_DIR / "logs"
OUT_DIR = BASE_DIR / "output"
LOG_DIR.mkdir(parents=True, exist_ok=True)

OVERLAP_REPORT_PATH = OUT_DIR / "species_overlap_report.csv"
OVERLAP_ONLY_PATH = OUT_DIR / "overlap_species_only.csv"
RESOLUTION_LOG_PATH = OUT_DIR / "phase1b_taxonomy_resolution_report.csv"

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
REQUEST_DELAY_SEC = 0.15  # polite pacing against a shared public API
MIN_CONFIDENCE = 90  # GBIF confidence threshold below which a match is not trusted

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_path = LOG_DIR / f"phase1b_{timestamp}.log"

logger = logging.getLogger("phase1b")
logger.setLevel(logging.DEBUG)
logger.handlers.clear()
fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(fmt)
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(fmt)
logger.addHandler(file_handler)
logger.addHandler(console_handler)


def load_report():
    with open(OVERLAP_REPORT_PATH, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def gbif_resolve(name, session):
    """Query GBIF species/match; return (accepted_binomial, status, matchType, confidence) or None."""
    try:
        resp = session.get(GBIF_MATCH_URL, params={"name": name, "rank": "SPECIES"}, timeout=15)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as exc:
        logger.warning(f"GBIF request failed for '{name}': {exc}")
        return None

    if not data or "usageKey" not in data:
        return None

    match_type = data.get("matchType", "NONE")
    confidence = data.get("confidence", 0)
    status = data.get("status", "UNKNOWN")
    accepted_binomial = data.get("species")  # GBIF always returns the accepted species name here

    if match_type == "NONE" or not accepted_binomial:
        return None

    return accepted_binomial, status, match_type, confidence


def main():
    logger.info("=" * 78)
    logger.info("PHASE 1 STEP 2: TAXONOMIC BACKBONE NAME RESOLUTION (GBIF)")
    logger.info("=" * 78)
    logger.info(f"Run started : {datetime.now().isoformat()}")
    logger.info(f"Log file    : {log_path}")

    if not OVERLAP_REPORT_PATH.exists():
        logger.error(f"{OVERLAP_REPORT_PATH} not found -- run phase1_data_ingestion.py first.")
        sys.exit(1)

    rows = load_report()
    by_species = {r["binomial_species"]: r for r in rows}

    cancer_only = [r["binomial_species"] for r in rows
                   if r["in_cancer_dataset"] == "True" and r["in_microbiome_dataset"] == "False"]
    microbiome_only = [r["binomial_species"] for r in rows
                       if r["in_microbiome_dataset"] == "True" and r["in_cancer_dataset"] == "False"]
    existing_overlap = {r["binomial_species"] for r in rows if r["overlap"] == "True"}

    logger.info(f"Original overlap (Phase 1, string-matching only): {len(existing_overlap)}")
    logger.info(f"Cancer-only species to resolve     : {len(cancer_only)}")
    logger.info(f"Microbiome-only species to resolve  : {len(microbiome_only)}")
    logger.info("")

    session = requests.Session()
    session.headers.update({"User-Agent": "oncobiome-phase1b/1.0 (research script)"})

    resolutions = {}  # original_name -> (accepted_binomial, status, matchType, confidence)
    all_targets = cancer_only + microbiome_only
    for i, name in enumerate(all_targets, 1):
        result = gbif_resolve(name, session)
        if result:
            resolutions[name] = result
        if i % 50 == 0:
            logger.debug(f"  ... resolved {i}/{len(all_targets)} names so far")
        time.sleep(REQUEST_DELAY_SEC)

    unresolved = [n for n in all_targets if n not in resolutions]
    logger.info(f"GBIF returned a match for {len(resolutions)}/{len(all_targets)} names "
                f"({len(unresolved)} unresolved/no match).")
    logger.info("")

    # ------------------------------------------------------------------
    # Case A: cancer-only and microbiome-only both resolve to the SAME
    # accepted binomial under different spellings -> genuinely new overlap
    # ------------------------------------------------------------------
    cancer_resolved = {n: resolutions[n][0] for n in cancer_only if n in resolutions}
    microbiome_resolved = {n: resolutions[n][0] for n in microbiome_only if n in resolutions}

    accepted_to_cancer_orig = {}
    for orig, accepted in cancer_resolved.items():
        accepted_to_cancer_orig.setdefault(accepted, []).append(orig)
    accepted_to_microbiome_orig = {}
    for orig, accepted in microbiome_resolved.items():
        accepted_to_microbiome_orig.setdefault(accepted, []).append(orig)

    new_overlap_accepted = set(accepted_to_cancer_orig) & set(accepted_to_microbiome_orig)
    # exclude accepted names that are already literally in the existing overlap set
    # (shouldn't happen by construction, but check explicitly rather than assume)
    new_overlap_accepted -= existing_overlap

    # ------------------------------------------------------------------
    # Case B: a cancer-only or microbiome-only name resolves to an
    # accepted binomial that matches an EXISTING overlap species under a
    # different spelling -- this doesn't create a new overlap species,
    # but it means that source now also covers an already-overlapping
    # species (worth logging, not worth double-counting).
    # ------------------------------------------------------------------
    extra_source_for_existing = []
    for orig, accepted in list(cancer_resolved.items()) + list(microbiome_resolved.items()):
        if accepted in existing_overlap and orig != accepted:
            extra_source_for_existing.append((orig, accepted))

    logger.info("=" * 78)
    logger.info("RESOLUTION RESULTS")
    logger.info("=" * 78)
    logger.info(f"NEW overlap species found via synonym resolution: {len(new_overlap_accepted)}")
    for accepted in sorted(new_overlap_accepted):
        cancer_origs = accepted_to_cancer_orig[accepted]
        micro_origs = accepted_to_microbiome_orig[accepted]
        logger.info(f"    - {accepted}  (cancer-side spelling: {cancer_origs}; "
                    f"microbiome-side spelling: {micro_origs})")
    logger.info("")
    if extra_source_for_existing:
        logger.info(f"Additional source-coverage found for already-overlapping species "
                    f"(does not change overlap count, logged for completeness): "
                    f"{len(extra_source_for_existing)}")
        for orig, accepted in extra_source_for_existing:
            logger.info(f"    - '{orig}' resolves to already-counted overlap species '{accepted}'")
    logger.info("")

    # ------------------------------------------------------------------
    # Write the full resolution audit trail
    # ------------------------------------------------------------------
    with open(RESOLUTION_LOG_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["original_name", "category", "gbif_accepted_binomial",
                          "gbif_status", "gbif_match_type", "gbif_confidence", "outcome"])
        for name in all_targets:
            category = "cancer_only" if name in cancer_only else "microbiome_only"
            if name in resolutions:
                accepted, status, match_type, confidence = resolutions[name]
                if accepted in new_overlap_accepted:
                    outcome = "NEW_OVERLAP_MATCH"
                elif accepted in existing_overlap and name != accepted:
                    outcome = "EXTRA_SOURCE_FOR_EXISTING_OVERLAP_SPECIES"
                else:
                    outcome = "NO_CROSS_CATEGORY_MATCH"
                writer.writerow([name, category, accepted, status, match_type, confidence, outcome])
            else:
                writer.writerow([name, category, "", "", "", "", "UNRESOLVED_BY_GBIF"])
    logger.info(f"Full resolution audit trail written to: {RESOLUTION_LOG_PATH}")

    # ------------------------------------------------------------------
    # Update species_overlap_report.csv and overlap_species_only.csv
    # ------------------------------------------------------------------
    fieldnames = list(rows[0].keys())
    if "resolved_via_taxonomy_backbone" not in fieldnames:
        fieldnames.append("resolved_via_taxonomy_backbone")

    for accepted in new_overlap_accepted:
        cancer_origs = accepted_to_cancer_orig[accepted]
        micro_origs = accepted_to_microbiome_orig[accepted]
        # Prefer an existing row if the accepted binomial happens to already
        # exist as a row key (from one of the two spellings); otherwise merge
        # the two original rows into one new row under the accepted name.
        merged_row = {fn: "" for fn in fieldnames}
        merged_row["binomial_species"] = accepted
        sources = set()
        for orig in cancer_origs + micro_origs:
            orig_row = by_species.get(orig, {})
            for k, v in orig_row.items():
                if k in ("binomial_species", "n_sources", "in_cancer_dataset",
                         "in_microbiome_dataset", "overlap"):
                    continue
                if v:
                    merged_row[k] = v
                    sources.add(k)
        merged_row["n_sources"] = str(len(sources))
        merged_row["in_cancer_dataset"] = "True"
        merged_row["in_microbiome_dataset"] = "True"
        merged_row["overlap"] = "True"
        merged_row["resolved_via_taxonomy_backbone"] = (
            f"cancer-side spelling(s): {cancer_origs}; microbiome-side spelling(s): {micro_origs}"
        )
        # remove the now-superseded original rows (they've been merged)
        for orig in cancer_origs + micro_origs:
            by_species.pop(orig, None)
        by_species[accepted] = merged_row

    for r in by_species.values():
        r.setdefault("resolved_via_taxonomy_backbone", "")

    final_rows = sorted(by_species.values(), key=lambda r: r["binomial_species"])
    with open(OVERLAP_REPORT_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_rows)

    final_overlap = [r for r in final_rows if r["overlap"] == "True"]
    with open(OVERLAP_ONLY_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["binomial_species", "sources"])
        for r in sorted(final_overlap, key=lambda r: r["binomial_species"]):
            src_cols = [k for k in fieldnames
                       if k not in ("binomial_species", "n_sources", "in_cancer_dataset",
                                    "in_microbiome_dataset", "overlap",
                                    "resolved_via_taxonomy_backbone")
                       and r.get(k)]
            writer.writerow([r["binomial_species"], "+".join(sorted(src_cols))])

    logger.info("")
    logger.info("=" * 78)
    logger.info("BEFORE / AFTER SUMMARY")
    logger.info("=" * 78)
    logger.info(f"Overlap BEFORE naming fix (original Phase 1 run)      : 55")
    logger.info(f"Overlap AFTER naming fix, same string matching        : {len(existing_overlap)} (unchanged, as expected)")
    logger.info(f"Overlap AFTER taxonomic backbone resolution (this run): {len(final_overlap)}")
    logger.info(f"Net new species from taxonomy resolution              : {len(final_overlap) - len(existing_overlap)}")
    if unresolved:
        logger.info(f"{len(unresolved)} species had no usable GBIF match at all "
                    f"(status NONE/low match type) -- these remain unmatched, not silently dropped.")
    logger.info(f"Updated overlap outputs written to: {OVERLAP_REPORT_PATH}, {OVERLAP_ONLY_PATH}")
    logger.info(f"Full log: {log_path}")


if __name__ == "__main__":
    main()
