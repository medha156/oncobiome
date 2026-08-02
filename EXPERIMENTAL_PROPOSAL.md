# Experimental Validation Proposal

## Status: design scoped, candidates not yet selected

This is a scoping document for a **follow-on study**, not something executed within this computational pipeline (`METHODOLOGY.md` Phase 7). **It is intentionally incomplete in one specific way**: the candidate species/metabolite claims this proposal would validate are supposed to come from Phase 6's flux-balance-analysis output (`output/phase6_fba_predictions.csv`). Phase 6 has not run — see `PIPELINE_EXECUTION_REPORT.md` for why (Phases 3-6 are blocked on missing raw sequencing data and missing bioinformatics tooling in this environment). Rather than invent plausible-sounding candidate species and metabolites to fill this section in, it is left as **TBD**, and this document should be revisited once Phase 5/6 produce real, logged findings.

---

## 1. Candidate species/metabolite claims

**TBD — pending real Phase 5 PGLS hits and Phase 6 FBA predictions.**

Once available, this section should name 1-3 candidates and state explicitly why each was prioritized over the rest of the Phase 6 hit list — per `METHODOLOGY.md` Phase 6/7, prioritize the clearest mechanistic story (a specific predicted metabolite with an established cancer-relevant mechanism: secondary bile acids, H2S, colibactin, or a butyrate deficit) over a merely statistically significant but mechanistically vague hit.

## 2. Proposed design

**Fecal microbiota transplant (FMT) into germ-free (gnotobiotic) mice**, comparing transplants from high-cancer-mortality-associated vs. low-cancer-mortality-associated source microbiomes (per `METHODOLOGY.md` Phase 7).

**Why germ-free FMT over the alternatives:**
- **vs. in vitro bioreactor systems (e.g. SHIME-type):** cheaper and avoid animal use, and are a reasonable, cheaper *pre-screen* for whether a donor community actually produces the Phase 6-predicted metabolites in a simplified gut-like environment — but they cannot capture host-immune or host-epithelial interactions, which is where the actual cancer-relevant phenotype (tumorigenesis, inflammation) plays out. **Recommended as an optional first step**, not a substitute for the mouse study.
- **vs. antibiotic-depletion mouse models:** cheaper and don't require gnotobiotic facilities, but leave a residual, incompletely-eliminated native community and residual antibiotic effects on host physiology — both of which confound interpretation of the transplanted phenotype. Germ-free mice avoid this confound entirely, which is why they're specified here despite the higher facility cost.

## 3. Pre-specified phenotype, endpoint, and sample size

To be filled in once a real candidate is selected, but written now with **Walter et al. 2020's critique of this exact experimental design explicitly in mind**: a systematic review of published human-microbiota-associated (HMA) rodent studies found 95% reported successful phenotype transfer — a rate the authors argue is implausibly high and likely reflects publication bias, small-effect-size over-interpretation, and confounding from incomplete/unstable colonization, not the gut microbiome being that reliably causal.

This proposal commits, in advance, to:
- Pre-registering the specific phenotype being tested (extending this project's Phase 0 preregistration practice to the animal experiment itself) before any transplant is performed.
- Using an adequately powered sample size per group, justified by an explicit power calculation at the design stage — not the field's typical underpowering that Walter et al. flag as part of why the 95% success rate is suspect.
- Reporting negative or partial-transfer results if that's what happens, not only positive ones.

## 4. Reporting commitment

This study will be reported using the **ARRIVE 2.0 guidelines** from the experimental design stage — sample size justification, randomization, and outcome reporting decided before data collection, not retrofitted at write-up time.

---

**This document is a design scaffold, not a finished proposal.** The parts that don't depend on which species/metabolite turns out to be the real Phase 6 hit (Sections 2-4) are complete. Section 1 needs real data before it can be filled in honestly.
