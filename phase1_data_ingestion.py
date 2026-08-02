#!/usr/bin/env python3
"""
Phase 1 -- data ingestion & species overlap check
====================================================
Project: Comparative analysis of mammalian gut microbiome functions
         associated with species-level cancer mortality

WHAT THIS SCRIPT DOES
----------------------
1. Loads (from a local cache in data/raw/) the species-level tables from the
   five source papers identified for this project.
2. Normalizes scientific names so they can be joined across datasets that
   used different naming conventions (underscores vs spaces, trinomial
   subspecies vs binomial species, etc).
3. Computes the overlap between species that have cancer-mortality data and
   species that have microbiome/metabolome data -- this overlap is the real
   usable sample size for every downstream step (PGLS, pathway analysis,
   FBA), so it needs to be known before any bioinformatics work starts.
4. Writes a full species-by-source matrix and an overlap-only species list
   to output/, and a detailed timestamped log to logs/.

DATA PROVENANCE (read before trusting numbers downstream)
-----------------------------------------------------------
- vincze2022_cancer_mortality.csv : REAL DATA. Downloaded verbatim from the
  Vincze et al. 2022 (Nature) authors' GitHub repo (data.csv), 191 species.
- manor2020_host_species.csv      : REAL DATA. Hand-extracted from Table 1
  of the Manor et al. 2020 (AEM) paper's full-text HTML, 77 species.
- compton2020_cancer_prevalence.csv,
  youngblut2020_host_species.csv,
  guo2022_host_species.csv         : REAL DATA. Extracted from supplementary
  Excel files the user downloaded manually and provided directly (see each
  source's 'manual_note' below), since pmc.ncbi.nlm.nih.gov, NCBI SRA, and
  MassIVE were not reachable from the sandbox network these were first
  attempted in.

USAGE
-----
    python phase1_data_ingestion.py

Re-runnable and idempotent: safe to run repeatedly as more source files
are added to data/raw/.
"""

import csv
import logging
import sys
from datetime import datetime
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
DATA_RAW = BASE_DIR / "data" / "raw"
LOG_DIR = BASE_DIR / "logs"
OUT_DIR = BASE_DIR / "output"

for d in (DATA_RAW, LOG_DIR, OUT_DIR):
    d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Logging -- console (INFO+) and a detailed timestamped file log (DEBUG+)
# --------------------------------------------------------------------------
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_path = LOG_DIR / f"phase1_{timestamp}.log"

logger = logging.getLogger("phase1")
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
# Data source registry
# --------------------------------------------------------------------------
DATA_SOURCES = {
    "vincze2022_cancer": {
        "citation": "Vincze et al. 2022, Nature 601:263-267 - Cancer risk across mammals",
        "source_url": "https://raw.githubusercontent.com/OrsolyaVincze/VinczeEtal2021Nature/main/data.csv",
        "local_file": DATA_RAW / "vincze2022_cancer_mortality.csv",
        "role": "cancer_mortality",
        "species_col": "Species",
        "expected_n": 191,
        "auto_fetchable": True,
    },
    "compton2020_cancer": {
        "citation": "Compton/Boddy et al. 2020, Evolution Medicine & Public Health - "
                    "Lifetime cancer prevalence and life history traits in mammals (PMC7652303)",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/instance/7652303/bin/eoaa015_supplementary_data.zip",
        "local_file": DATA_RAW / "compton2020_cancer_prevalence.csv",
        "role": "cancer_mortality",
        "species_col": "species_name",
        "expected_n": 37,
        "auto_fetchable": True,
        "manual_note": (
            "RESOLVED 2026-08-01: user manually downloaded 'supplementary_tables_V2.xlsx' "
            "(sheet S2_LHtable) and provided it directly, since pmc.ncbi.nlm.nih.gov is not "
            "reachable from this environment's network allow-list. Extracted to local_file."
        ),
    },
    "manor2020_microbiome": {
        "citation": "Manor et al. 2020, Applied and Environmental Microbiology 86:e01864-20 - "
                    "Multi-omics Approaches To Decipher the Impact of Diet and Host Physiology "
                    "on the Mammalian Gut Microbiome",
        "source_url": "https://journals.asm.org/doi/10.1128/aem.01864-20",
        "local_file": DATA_RAW / "manor2020_host_species.csv",
        "role": "microbiome",
        "species_col": "species",
        "expected_n": 77,
        "auto_fetchable": True,
        "sra_bioprojects": ["PRJNA545289", "PRJNA545214"],
    },
    "youngblut2020_microbiome": {
        "citation": "Youngblut et al. 2020, mSystems 5:e01045-20 - Large-Scale Metagenome "
                    "Assembly Reveals Novel Animal-Associated Microbial Genomes...",
        "source_url": "https://github.com/leylabmpi/animal_gut_metagenome_assembly",
        "local_file": DATA_RAW / "youngblut2020_host_species.csv",
        "role": "microbiome",
        "species_col": "species",
        "expected_n": 180,
        "auto_fetchable": True,
        "manual_note": (
            "RESOLVED 2026-08-01: user manually downloaded the paper's Table S1 "
            "(mSystems.01045-20-st001.xlsx, sheet 'Table S1A') and provided it directly. "
            "NOTE: this 180-species count spans ALL 5 host classes sampled (Mammalia, "
            "Aves, Actinopterygii, Reptilia, Amphibia), not mammals only -- only 103 of "
            "the 180 are Mammalia, which is the subset relevant to the cancer-mortality "
            "overlap. It also does NOT include species covered only by the 14 external "
            "NCBI BioProjects listed in Table S1B (pig, dog, macaque, minke whale, cod, "
            "house cat, goose, woodrat/grouse, kakapo, chicken, pangolin, capuchin, black "
            "rhinoceros) since those are given as common names only, not resolved to "
            "binomial species names -- a manual name resolution follow-up."
        ),
    },
    "guo2022_metabolome": {
        "citation": "Gregor/Guo et al. 2022, The ISME Journal 16:1262-1274 - Mammalian gut "
                    "metabolomes mirror microbiome composition and host phylogeny",
        "source_url": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000086131",
        "local_file": DATA_RAW / "guo2022_host_species.csv",
        "role": "microbiome",
        "species_col": "Taxonomy",
        "expected_n": 25,
        "auto_fetchable": True,
        "manual_note": (
            "RESOLVED 2026-08-01: user manually downloaded the paper's Supplementary Data "
            "2 (41396_2021_1152_moesm2_esm.xlsx, sheet 'Feces Sampling Table') and provided "
            "it directly. Some entries are subspecies trinomials (e.g. 'Canis lupus arabs', "
            "'Varecia variegata rubra') which the normalizer below collapses to binomial."
        ),
    },
}


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_source(key, cfg):
    logger.info(f"--- Loading source: {key} ---")
    logger.info(f"Citation   : {cfg['citation']}")
    logger.info(f"Source URL : {cfg['source_url']}")
    logger.debug(f"Local file : {cfg['local_file']}")

    if not cfg["local_file"].exists():
        if cfg["auto_fetchable"]:
            logger.error(
                f"Expected local cache at {cfg['local_file']} but it is missing. "
                f"This source is normally auto-fetchable -- re-fetch from "
                f"{cfg['source_url']} and save it there before continuing."
            )
        else:
            logger.warning(f"No local file found for '{key}'. Marking UNAVAILABLE this run.")
            logger.warning(f"Manual retrieval needed: {cfg.get('manual_note', 'n/a')}")
        return None

    species = set()
    rows_read = 0
    blank_rows = 0
    with open(cfg["local_file"], newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if cfg["species_col"] not in (reader.fieldnames or []):
            logger.error(
                f"Column '{cfg['species_col']}' not found in {cfg['local_file']}. "
                f"Columns present: {reader.fieldnames}"
            )
            return None
        for row in reader:
            rows_read += 1
            raw_name = (row.get(cfg["species_col"]) or "").strip()
            if not raw_name:
                blank_rows += 1
                continue
            species.add(raw_name)

    logger.info(f"Rows read            : {rows_read}")
    if blank_rows:
        logger.warning(f"Blank species values : {blank_rows} row(s) skipped")
    logger.info(f"Unique species names  : {len(species)}")

    expected_n = cfg.get("expected_n")
    if expected_n and len(species) != expected_n:
        logger.warning(
            f"Expected ~{expected_n} species per the source paper, got {len(species)}. "
            f"Check for duplicate or malformed rows before trusting downstream counts."
        )
    elif expected_n:
        logger.info(f"Species count matches paper's reported n ({expected_n}). OK.")

    return species


# --------------------------------------------------------------------------
# Name normalization
# --------------------------------------------------------------------------
def normalize_species_name(name):
    """
    Standardize a scientific name for cross-dataset joining.

    - underscores -> spaces, whitespace collapsed
    - Genus capitalized, species epithet lowercased
    - Trinomials (subspecies, e.g. 'Equus africanus asinus') are collapsed
      to their binomial ('Equus africanus') for matching purposes, since
      the cancer-mortality datasets are species-level, not subspecies-level.
      The full original name is preserved in the per-source columns of the
      output report, so no information is discarded -- only the join key
      is coarsened.

    KNOWN LIMITATION: this is string-based normalization only. It does NOT
    resolve taxonomic synonyms (e.g. a species renamed between genera in a
    revision). A production version should cross-check names against a
    taxonomy backbone (NCBI Taxonomy / GBIF), which requires network access
    this environment did not have. Treat near-miss species (present in one
    dataset under a slightly different name) as a manual follow-up.
    """
    cleaned = name.replace("_", " ").strip()
    cleaned = " ".join(cleaned.split())
    parts = cleaned.split(" ")
    if len(parts) >= 2 and parts[0] and parts[1]:
        binomial = f"{parts[0][0].upper()}{parts[0][1:].lower()} {parts[1].lower()}"
    else:
        binomial = cleaned
    return cleaned, binomial


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    logger.info("=" * 78)
    logger.info("PHASE 1: DATA INGESTION & SPECIES OVERLAP CHECK")
    logger.info("=" * 78)
    logger.info(f"Run started : {datetime.now().isoformat()}")
    logger.info(f"Log file    : {log_path}")
    logger.info(f"Data sources: {len(DATA_SOURCES)} registered")
    logger.info("")

    raw_species_by_source = {}
    for key, cfg in DATA_SOURCES.items():
        raw_species_by_source[key] = load_source(key, cfg)
        logger.info("")

    logger.info("=" * 78)
    logger.info("SOURCE AVAILABILITY SUMMARY")
    logger.info("=" * 78)
    available = [k for k, v in raw_species_by_source.items() if v is not None]
    unavailable = [k for k, v in raw_species_by_source.items() if v is None]
    logger.info(f"Loaded this run ({len(available)}/{len(DATA_SOURCES)}): {available}")
    logger.info(f"Needs manual download ({len(unavailable)}/{len(DATA_SOURCES)}): {unavailable}")
    logger.info("")

    # ---------------------------------------------------------- normalize --
    logger.info("=" * 78)
    logger.info("NORMALIZING SPECIES NAMES")
    logger.info("=" * 78)
    normalized = {}  # binomial -> {source_key: raw_name}
    for key, species_set in raw_species_by_source.items():
        if species_set is None:
            continue
        n_before = len(species_set)
        collisions = 0
        for raw in species_set:
            _, binomial = normalize_species_name(raw)
            normalized.setdefault(binomial, {})
            if key in normalized[binomial]:
                collisions += 1
            normalized[binomial][key] = raw
        n_after = len({b for b, srcs in normalized.items() if key in srcs})
        logger.info(
            f"{key}: {n_before} raw name(s) -> {n_after} unique binomial(s) "
            f"({collisions} subspecies collapsed onto an existing binomial)"
        )
    logger.info(f"Total unique binomial species names across ALL loaded sources: {len(normalized)}")
    logger.info("")

    # ---------------------------------------------------------- overlap ----
    cancer_keys = [k for k, v in DATA_SOURCES.items() if v["role"] == "cancer_mortality"]
    microbiome_keys = [k for k, v in DATA_SOURCES.items() if v["role"] == "microbiome"]

    cancer_species = {b for b, srcs in normalized.items() if any(k in srcs for k in cancer_keys)}
    microbiome_species = {b for b, srcs in normalized.items() if any(k in srcs for k in microbiome_keys)}
    overlap_species = cancer_species & microbiome_species

    logger.info("=" * 78)
    logger.info("OVERLAP CHECK: cancer-mortality data x microbiome/metabolome data")
    logger.info("=" * 78)
    loaded_cancer = [k for k in cancer_keys if raw_species_by_source[k] is not None]
    loaded_microbiome = [k for k in microbiome_keys if raw_species_by_source[k] is not None]
    logger.info(f"Species with cancer-mortality data     : {len(cancer_species)}  (from: {loaded_cancer})")
    logger.info(f"Species with microbiome/metabolome data: {len(microbiome_species)}  (from: {loaded_microbiome})")
    logger.info(f"*** SPECIES IN BOTH (usable for downstream PGLS / pathway / FBA work): {len(overlap_species)} ***")
    for sp in sorted(overlap_species):
        srcs = normalized[sp]
        src_tags = "+".join(sorted(srcs.keys()))
        logger.info(f"    - {sp:<28} [{src_tags}]")
    logger.info("")

    if unavailable:
        logger.warning(
            f"This overlap count is a LOWER BOUND: {len(unavailable)} source(s) could not "
            f"be auto-loaded this run ({unavailable}) and will likely add more overlapping "
            f"species once retrieved manually -- see per-source manual_note above."
        )
        logger.info("")

    # ---------------------------------------------------------- outputs ----
    report_path = OUT_DIR / "species_overlap_report.csv"
    fieldnames = ["binomial_species", "n_sources"] + list(DATA_SOURCES.keys()) + [
        "in_cancer_dataset",
        "in_microbiome_dataset",
        "overlap",
    ]
    with open(report_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for binomial in sorted(normalized.keys()):
            srcs = normalized[binomial]
            row = {"binomial_species": binomial, "n_sources": len(srcs)}
            for key in DATA_SOURCES.keys():
                row[key] = srcs.get(key, "")
            row["in_cancer_dataset"] = binomial in cancer_species
            row["in_microbiome_dataset"] = binomial in microbiome_species
            row["overlap"] = binomial in overlap_species
            writer.writerow(row)
    logger.info(f"Full species-by-source matrix written to : {report_path}")

    overlap_path = OUT_DIR / "overlap_species_only.csv"
    with open(overlap_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["binomial_species", "sources"])
        for sp in sorted(overlap_species):
            writer.writerow([sp, "+".join(sorted(normalized[sp].keys()))])
    logger.info(f"Overlap-only species list written to     : {overlap_path}")

    logger.info("")
    logger.info("=" * 78)
    logger.info("PHASE 1 COMPLETE")
    logger.info("=" * 78)
    logger.info(f"Full detailed log saved to: {log_path}")
    if unavailable:
        logger.info(f"Next step: manually retrieve {len(unavailable)} unavailable source(s) "
                    "(see warnings above), drop them in data/raw/, and re-run to widen the overlap.")
    else:
        logger.info("All registered sources loaded. Next step: taxonomic name resolution "
                    "(e.g. GBIF/NCBI Taxonomy) to catch near-miss synonyms this string-based "
                    "normalizer can't, then move to Phase 2 (bioinformatics processing).")


if __name__ == "__main__":
    main()
