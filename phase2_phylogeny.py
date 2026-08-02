#!/usr/bin/env python3
"""
Phase 2 -- species-level time-calibrated phylogeny
====================================================
Project: Comparative analysis of mammalian gut microbiome functions
         associated with species-level cancer mortality

WHAT THIS SCRIPT DOES
----------------------
1. Reads the cancer x microbiome overlap list produced by Phase 1
   (output/overlap_species_only.csv) -- 57 species as of the taxonomy
   backbone resolution pass (phase1b_taxonomy_resolution.py).
2. Submits that species list to TimeTree (timetree.org) -- the same
   resource Youngblut et al. 2020 used -- and retrieves a pruned,
   time-calibrated Newick tree spanning exactly those species.
3. Cross-checks every input species against the resulting tree tips.
   This is deliberately done now, before any PGLS/statistics, because
   TimeTree's own name resolution surfaces exactly the kind of
   cross-dataset taxonomic naming problems Phase 1's normalizer cannot
   catch (synonyms, genus reassignments, insufficient data for a taxon).
4. Writes the raw tree, a per-species match report, and a timestamped
   log.

TIMETREE VS. VERTLIFE/UPHAM -- WHY TIMETREE WAS USED
------------------------------------------------------
METHODOLOGY.md Phase 2 names the VertLife/Upham et al. 2019 mammal
supertree as the *preferred* source (it ships as a credible set of
trees, not a single point estimate, which lets PGLS be re-run across
tree replicates to check sensitivity to phylogenetic uncertainty).
This was checked for feasibility before defaulting to TimeTree:

  - The full Dryad data package (doi:10.5061/dryad.tb03d03) is 5.5 GB --
    impractical to fetch in bulk just to extract one species subset.
  - The smaller, specifically relevant file (Data_S4_patchClade_results_
    and_MCC.zip, ~4.5 MB, likely containing the maximum-clade-credibility
    single tree) could not actually be downloaded from this environment:
    the Dryad REST API's file-download endpoint requires an OAuth bearer
    token (401 Unauthorized), and the public web "file_stream" download
    link is behind a Cloudflare JS challenge ("Validating..." page) that
    a plain HTTP client cannot pass.
  - The interactive VertLife species-subsetting tool (vertlife.org/
    phylosubsets) does load, but is not a simple stateless request/
    response flow like TimeTree's -- it did not yield a quick,
    reasonably-scoped path to a scripted subset download either.

Given this, TimeTree was used as the documented fallback, per
METHODOLOGY.md's own instruction to do so explicitly rather than
silently substitute. **This is a real limitation, not a stylistic
choice**: it means this project's PGLS step (Phase 5) will use a single
point-estimate tree and cannot re-run models across a credible set of
trees to quantify sensitivity to phylogenetic uncertainty, the way an
Upham-et-al.-based analysis could.

DATA PROVENANCE / IMPORTANT CAVEAT
-----------------------------------
TimeTree (timetree.org) has no documented public REST API. The "Load a
List" -> "Prune Tree" feature on the website is implemented as a pair of
undocumented AJAX endpoints used by their own front-end JavaScript:

    POST /ajax/prune/load_names/         (multipart file upload, field "file")
    POST /ajax/newick/prunetree/download  (form field export=newick)

This script drives those two endpoints directly with a `requests`
session (reverse-engineered from the site's public/js/app.js -- no
authentication or private access involved, same as what a browser does
when a user manually uses the "Load a List" tool). This is NOT a stable,
versioned API: TimeTree can change these endpoints at any time without
notice, which would break this script. Treat it the same way Phase 1
treats its manually-downloaded sources -- re-verify if it stops working,
and if TimeTree changes their site, the fallback is to use the
interactive "Load a List" web tool by hand and drop the resulting
Newick file at RAW_TREE_PATH.

USAGE
-----
    python phase2_phylogeny.py

Re-runnable and idempotent: safe to re-run if the Phase 1 overlap list
changes (e.g. once more species are added / names reconciled).
"""

import csv
import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import requests

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
DATA_RAW = BASE_DIR / "data" / "raw"
DATA_PROCESSED = BASE_DIR / "data" / "processed"
PHYLO_DIR = DATA_RAW / "phylogeny"
LOG_DIR = BASE_DIR / "logs"
OUT_DIR = BASE_DIR / "output"

for d in (PHYLO_DIR, DATA_PROCESSED, LOG_DIR, OUT_DIR):
    d.mkdir(parents=True, exist_ok=True)

OVERLAP_LIST_PATH = OUT_DIR / "overlap_species_only.csv"
SUBMITTED_LIST_PATH = PHYLO_DIR / "timetree_submitted_species_list.txt"
RAW_TREE_PATH = PHYLO_DIR / "timetree_prunetree_raw.nwk"
RAW_HTML_PATH = PHYLO_DIR / "timetree_prunetree_response.html"
MATCH_REPORT_PATH = OUT_DIR / "phylogeny_species_match_report.csv"
PGLS_TREE_PATH = PHYLO_DIR / "timetree_pgls_ready.nwk"
PGLS_SPECIES_LIST_PATH = OUT_DIR / "pgls_species_list.csv"
# Canonical path later pipeline phases (Phase 5 prompt) expect the final
# pruned, PGLS-ready tree at -- identical content to PGLS_TREE_PATH, just
# the stable name/location downstream phases read from.
FINAL_TREE_PATH = DATA_PROCESSED / "overlap_species_tree.nwk"

TIMETREE_BASE = "http://www.timetree.org"
LOAD_NAMES_URL = f"{TIMETREE_BASE}/ajax/prune/load_names/"
DOWNLOAD_NEWICK_URL = f"{TIMETREE_BASE}/ajax/newick/prunetree/download"

# --------------------------------------------------------------------------
# Logging -- console (INFO+) and a detailed timestamped file log (DEBUG+)
# --------------------------------------------------------------------------
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_path = LOG_DIR / f"phase2_{timestamp}.log"

logger = logging.getLogger("phase2")
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


# --------------------------------------------------------------------------
# Name normalization (mirrors phase1_data_ingestion.normalize_species_name)
# --------------------------------------------------------------------------
def normalize_binomial(name):
    cleaned = name.replace("_", " ").strip()
    cleaned = " ".join(cleaned.split())
    parts = cleaned.split(" ")
    if len(parts) >= 2 and parts[0] and parts[1]:
        return f"{parts[0][0].upper()}{parts[0][1:].lower()} {parts[1].lower()}"
    return cleaned


def load_overlap_species():
    species = []
    with open(OVERLAP_LIST_PATH, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            species.append(row["binomial_species"].strip())
    return species


# --------------------------------------------------------------------------
# TimeTree fetch
# --------------------------------------------------------------------------
def fetch_timetree_prunetree(species_list):
    session = requests.Session()
    session.headers.update({"User-Agent": "oncobiome-phase2/1.0 (research script)"})

    logger.info(f"Establishing session with {TIMETREE_BASE} ...")
    resp = session.get(TIMETREE_BASE, timeout=30)
    resp.raise_for_status()
    logger.debug(f"Homepage GET status: {resp.status_code}")

    species_blob = "\n".join(species_list) + "\n"
    SUBMITTED_LIST_PATH.write_text(species_blob, encoding="utf-8")
    logger.info(f"Submitted species list saved to: {SUBMITTED_LIST_PATH}")

    logger.info(f"POST {len(species_list)} species names to {LOAD_NAMES_URL} ...")
    files = {"file": ("species_list.txt", species_blob, "text/plain")}
    resp = session.post(LOAD_NAMES_URL, files=files, timeout=60)
    resp.raise_for_status()
    RAW_HTML_PATH.write_text(resp.text, encoding="utf-8")
    logger.debug(f"Prune-tree widget HTML saved to: {RAW_HTML_PATH}")

    if "search-error" in resp.text or "No Results" in resp.text:
        logger.error(
            "TimeTree returned an error page instead of a tree widget. "
            f"See {RAW_HTML_PATH} for the raw response. Aborting fetch."
        )
        return None, resp.text

    logger.info(f"POST export=newick to {DOWNLOAD_NEWICK_URL} ...")
    resp = session.post(DOWNLOAD_NEWICK_URL, data={"export": "newick"}, timeout=60)
    resp.raise_for_status()
    newick = resp.text.strip()

    if not newick or "(" not in newick:
        logger.error(
            "Newick download did not return a parseable tree "
            f"(got {len(newick)} bytes). Aborting."
        )
        return None, resp.text

    RAW_TREE_PATH.write_text(newick + "\n", encoding="utf-8")
    logger.info(f"Raw Newick tree saved to: {RAW_TREE_PATH} ({len(newick)} bytes)")

    return newick, RAW_HTML_PATH.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------
def extract_tip_labels(newick):
    """Pull leaf labels (Genus_species, underscore form) out of a Newick string."""
    tips = re.findall(r"[,(]([A-Za-z][A-Za-z_.\-]+):", newick)
    return sorted(set(tips))


def extract_unresolved_notes(html):
    """
    Parse the '#unresolved-names' box TimeTree renders when it had to
    substitute or drop a submitted name, e.g.:
        "Gazella subgutturosa (replaced with Gazella dorcas)"
        "Felis silvestris (insufficient data in TimeTree to place this taxon)"
    Returns {input_binomial: note_string}.
    """
    notes = {}
    match = re.search(
        r'id="unresolved-names">(.*?)</div>', html, flags=re.DOTALL
    )
    if not match:
        return notes
    block = match.group(1)
    for line in block.split("<br/>"):
        line = line.strip()
        if not line:
            continue
        m = re.match(r"([A-Za-z]+ [a-z][a-z\-]*)\s*\((.+)\)$", line)
        if m:
            notes[normalize_binomial(m.group(1))] = m.group(2).strip()
    return notes


def extract_leaf_count(html):
    match = re.search(r'id="prunetree-leaf-count">\s*(\d+)', html, flags=re.DOTALL)
    return int(match.group(1)) if match else None


# --------------------------------------------------------------------------
# Cross-check tree tips against the input overlap list
# --------------------------------------------------------------------------
def load_overlap_sources():
    """binomial_species -> sources string, from Phase 1's overlap output."""
    sources = {}
    with open(OVERLAP_LIST_PATH, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sources[row["binomial_species"].strip()] = row["sources"].strip()
    return sources


def finalize_pgls_tree(newick, rows):
    """
    Produce a PGLS-ready tree + species list from the raw TimeTree pull and
    its match report:
      - silent_genus_rename tips are relabeled back to the binomial used in
        the Phase 1 data tables (join-key fix only; topology/branch lengths
        untouched), so PGLS tools can match tree tips to trait-table rows
        by exact string.
      - dropped_no_data / missing_unexplained species are excluded from the
        final species list (nothing to rename -- they were never placed on
        the tree in the first place).
    """
    rename_map = {}
    for r in rows:
        if r["status"] == "silent_genus_rename":
            target = normalize_binomial(r["input_species"]).replace(" ", "_")
            rename_map[r["tree_tip"]] = target

    pgls_newick = newick
    for old_tip, new_tip in rename_map.items():
        pgls_newick = re.sub(
            rf"(?<=[,(]){re.escape(old_tip)}(?=:)", new_tip, pgls_newick
        )
    PGLS_TREE_PATH.write_text(pgls_newick + "\n", encoding="utf-8")
    logger.info(f"PGLS-ready tree (tips relabeled) written to: {PGLS_TREE_PATH}")
    for old_tip, new_tip in rename_map.items():
        logger.info(f"    relabeled tip: {old_tip} -> {new_tip}")
    FINAL_TREE_PATH.write_text(pgls_newick + "\n", encoding="utf-8")
    logger.info(f"Same tree also written to canonical path: {FINAL_TREE_PATH} "
                f"(this is what Phase 5 reads)")

    sources = load_overlap_sources()
    excluded = [r for r in rows if r["status"] in ("dropped_no_data", "missing_unexplained")]
    final_rows = []
    for r in rows:
        if r["status"] in ("dropped_no_data", "missing_unexplained"):
            continue
        binomial = normalize_binomial(r["input_species"])
        note = r["detail"] if r["status"] != "exact_match" else r["detail"]
        final_rows.append(
            {
                "binomial_species": binomial,
                "sources": sources.get(r["input_species"], ""),
                "note": note,
            }
        )

    with open(PGLS_SPECIES_LIST_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["binomial_species", "sources", "note"])
        writer.writeheader()
        writer.writerows(final_rows)
    logger.info(
        f"Final PGLS species list ({len(final_rows)} species) written to: "
        f"{PGLS_SPECIES_LIST_PATH}"
    )
    if excluded:
        logger.info(
            f"Excluded from PGLS sample ({len(excluded)}): "
            f"{[r['input_species'] for r in excluded]}"
        )
    return final_rows


def build_match_report(input_species, tip_labels, unresolved_notes):
    tip_binomials = {normalize_binomial(t): t for t in tip_labels}
    # epithet -> list of tip binomials sharing that species epithet, for
    # detecting silent genus reassignments TimeTree doesn't flag explicitly.
    epithet_index = {}
    for binomial in tip_binomials:
        genus, _, epithet = binomial.partition(" ")
        epithet_index.setdefault(epithet, []).append(binomial)

    rows = []
    for sp in input_species:
        norm = normalize_binomial(sp)
        note = unresolved_notes.get(norm, "")

        if norm in tip_binomials:
            status = "exact_match"
            tree_tip = tip_binomials[norm]
            detail = note or ""
        else:
            genus, _, epithet = norm.partition(" ")
            candidates = [b for b in epithet_index.get(epithet, []) if b != norm]
            if candidates:
                status = "silent_genus_rename"
                tree_tip = tip_binomials[candidates[0]]
                detail = (
                    f"Not found as '{norm}'; matched by species epithet to "
                    f"tree tip '{tree_tip}' (TimeTree did not flag this in its "
                    f"own unresolved-names list -- likely a taxonomic revision "
                    f"not reflected in the source dataset's naming)."
                )
            elif note:
                status = "dropped_no_data"
                tree_tip = ""
                detail = note
            else:
                status = "missing_unexplained"
                tree_tip = ""
                detail = "Not present in tree and not mentioned in TimeTree's unresolved-names list."

        rows.append(
            {
                "input_species": sp,
                "tree_tip": tree_tip,
                "status": status,
                "detail": detail,
            }
        )
    return rows


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    logger.info("=" * 78)
    logger.info("PHASE 2: SPECIES-LEVEL TIME-CALIBRATED PHYLOGENY (TimeTree)")
    logger.info("=" * 78)
    logger.info(f"Run started : {datetime.now().isoformat()}")
    logger.info(f"Log file    : {log_path}")
    logger.info("")
    logger.info("Tree source decision: VertLife/Upham et al. 2019 mammal supertree was the "
                "preferred source (credible tree set, not a single point estimate -- see "
                "METHODOLOGY.md Phase 2). Checked for feasibility: full Dryad package is "
                "5.5 GB; the smaller MCC-tree file is blocked by a Cloudflare JS challenge on "
                "direct download (and the Dryad REST API requires an OAuth bearer token this "
                "environment doesn't have). Falling back to TimeTree, per METHODOLOGY.md's "
                "explicit instruction to note the fallback and why. LIMITATION: this means "
                "Phase 5 will use a single point-estimate tree, not a credible tree-set, and "
                "cannot quantify PGLS sensitivity to phylogenetic uncertainty the way an "
                "Upham-et-al.-based analysis could.")
    logger.info("")

    if not OVERLAP_LIST_PATH.exists():
        logger.error(
            f"{OVERLAP_LIST_PATH} not found -- run phase1_data_ingestion.py first."
        )
        sys.exit(1)

    input_species = load_overlap_species()
    logger.info(f"Loaded {len(input_species)} overlap species from Phase 1 output.")
    logger.info("")

    try:
        newick, html = fetch_timetree_prunetree(input_species)
    except requests.RequestException as exc:
        logger.error(f"Network error contacting TimeTree: {exc}")
        sys.exit(1)

    if newick is None:
        logger.error("Could not obtain a tree from TimeTree this run. See log above.")
        sys.exit(1)

    tip_labels = extract_tip_labels(newick)
    unresolved_notes = extract_unresolved_notes(html)
    leaf_count = extract_leaf_count(html)

    logger.info("")
    logger.info("=" * 78)
    logger.info("TREE SUMMARY")
    logger.info("=" * 78)
    logger.info(f"Species submitted        : {len(input_species)}")
    logger.info(f"Leaf nodes in tree       : {leaf_count if leaf_count is not None else len(tip_labels)}")
    if unresolved_notes:
        logger.info(f"TimeTree-flagged name issues ({len(unresolved_notes)}):")
        for sp, note in unresolved_notes.items():
            logger.info(f"    - {sp}: {note}")
    logger.info("")

    rows = build_match_report(input_species, tip_labels, unresolved_notes)

    status_counts = {}
    for r in rows:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1

    logger.info("=" * 78)
    logger.info("SPECIES <-> TREE TIP CROSS-CHECK")
    logger.info("=" * 78)
    for status, count in sorted(status_counts.items()):
        logger.info(f"{status:<22}: {count}")
    logger.info("")
    for r in rows:
        if r["status"] != "exact_match":
            logger.warning(f"[{r['status']}] {r['input_species']} -> {r['detail']}")
    logger.info("")

    with open(MATCH_REPORT_PATH, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["input_species", "tree_tip", "status", "detail"]
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Per-species match report written to: {MATCH_REPORT_PATH}")
    logger.info("")

    logger.info("=" * 78)
    logger.info("FINALIZING PGLS-READY TREE + SPECIES LIST")
    logger.info("=" * 78)
    finalize_pgls_tree(newick, rows)

    usable = sum(1 for r in rows if r["status"] in ("exact_match", "silent_genus_rename"))
    logger.info("")
    logger.info("=" * 78)
    logger.info("PHASE 2 COMPLETE")
    logger.info("=" * 78)
    logger.info(
        f"{usable}/{len(input_species)} overlap species are placed on the tree "
        f"and usable for downstream PGLS."
    )
    dropped = [r["input_species"] for r in rows if r["status"] == "dropped_no_data"]
    if dropped:
        logger.info(
            f"{len(dropped)} species have NO tree placement and must be excluded "
            f"from PGLS (or manually grafted onto the tree at a conservative node): "
            f"{dropped}"
        )
    renamed = [r["input_species"] for r in rows if r["status"] == "silent_genus_rename"]
    if renamed:
        logger.info(
            f"{len(renamed)} species matched only after allowing for a genus-level "
            f"taxonomic revision -- flag these for the taxonomy cross-check step: {renamed}"
        )
    logger.info(f"Full detailed log saved to: {log_path}")
    logger.info(
        "Next step: use RAW_TREE_PATH as the PGLS phylogeny, substituting tip "
        "labels per the match report so they align with the species names used "
        "in the cancer-mortality and microbiome tables."
    )


if __name__ == "__main__":
    main()
