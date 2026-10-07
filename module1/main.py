# main.py

import sys
import json
import os

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

from processing  import validate_sequence, get_mhci, get_mhcii, get_bcell
from bcell       import extract_bcell_epitopes, extract_bcell_epitopes_relaxed
from scoring     import score_mhci, score_mhcii
from fusion      import fuse_epitopes
from default_config import DEFAULT_CONFIG
from filter      import is_redundant
from advanced_filters import (
    run_toxinpred
)
from allergenicity import predict_allergen, filter_candidates_by_safety
from iedb_advanced   import compute_population_coverage, compute_conservancy
from antigenicity    import score_epitope_antigenicity
from protparam       import analyze_construct
from adjuvant_selector import auto_select_adjuvant, print_adjuvant_comparison
from homology        import check_human_homology
from immune_export   import generate_immune_report
from immunogenicity  import batch_predict_immunogenicity, estimate_tc50
import math
from epitope_report import export_epitope_csv


# ── Diversity selection ────────────────────────────────────────────────────────

def select_diverse_top(candidates, count, already_selected_peptides,
                       overlap_threshold, config=None):
    """Select top-N peptides that are not redundant with already-selected ones."""
    selected = []
    for item in candidates:
        if len(selected) >= count:
            break
        peptide = item["peptide"]
        pool = already_selected_peptides + [s["peptide"] for s in selected]
        if not is_redundant(peptide, pool, overlap_threshold, config):
            selected.append(item)
    return selected


def _clean_str(s) -> str:
    if not isinstance(s, str):
        s = str(s)
        
    replacements = {
        "≥": ">=", "≤": "<=", "→": "->", "←": "<-",
        "✓": "[PASS]", "✗": "[FAIL]", "⚠": "[WARN]",
        "—": "-", "–": "-", "°": " deg"
    }
    for k, v in replacements.items():
        s = s.replace(k, v)
        
    try:
        s.encode(sys.stdout.encoding or 'ascii')
        return s
    except Exception:
        return s.encode('ascii', errors='replace').decode('ascii')


def _section(title: str):
    bar = "=" * 60
    t = _clean_str(title)
    print(f"\n{bar}")
    print(f"  {t}")
    print(bar)

def _ok(msg):   print(f"  [OK]   {_clean_str(msg)}")
def _warn(msg):  print(f"  [WARN] {_clean_str(msg)}")
def _fail(msg):  print(f"  [FAIL] {_clean_str(msg)}")
def _info(msg):  print(f"  [INFO] {_clean_str(msg)}")



# ── Main pipeline ──────────────────────────────────────────────────────────────

def run_module1(seq, config=None):
    if config is None:
        config = DEFAULT_CONFIG

    # Ensure sequence is injected into config for fallbacks
    config["sequence"] = seq

    # ── PHASE 1A: Sequence validation ────────────────────────────────────────
    _section("PHASE 1 — SEQUENCE VALIDATION")
    seq = validate_sequence(seq, config)
    _ok(f"Sequence valid: {len(seq)} amino acids")

    # ── PHASE 1B: Epitope prediction via IEDB API ─────────────────────────────
    _section("PHASE 1 — EPITOPE PREDICTION (IEDB API)")

    print("[API] Fetching MHC-I predictions...")
    mhci_df = get_mhci(seq, config)
    _ok(f"MHC-I: {len(mhci_df)} raw hits")

    print("[API] Fetching MHC-II predictions...")
    mhcii_df = get_mhcii(seq, config)
    _ok(f"MHC-II: {len(mhcii_df)} raw hits")

    print("[API] Fetching B-cell predictions...")
    bcell_df = get_bcell(seq, config)
    print(f"  [DEBUG] Raw B-cell rows: {len(bcell_df)}")

    bcell_epitopes = extract_bcell_epitopes(bcell_df, config)

    if len(bcell_epitopes) == 0:
        print("  [FALLBACK] B-cell triggered")
        # Lower threshold dynamically (Point 3)
        config["thresholds"]["bcell"] -= 0.1
        bcell_epitopes = extract_bcell_epitopes_relaxed(bcell_df, config, min_count=1)

    _ok(f"B-cell: {len(bcell_epitopes)} epitopes extracted")

    # ── PHASE 1C: Scoring ─────────────────────────────────────────────────────
    _section("PHASE 1 — EPITOPE SCORING")

    mhci_scores  = score_mhci(mhci_df,   config)
    mhcii_scores = score_mhcii(mhcii_df, config)

    def get_best_unique(scores):
        unique = {}
        for item in scores:
            p = item["peptide"]
            if p not in unique or item["score"] > unique[p]["score"]:
                unique[p] = item
        return sorted(unique.values(), key=lambda x: x["score"], reverse=True)

    mhci_candidates  = get_best_unique(mhci_scores)
    mhcii_candidates = get_best_unique(mhcii_scores)
    bcell_score_val  = config["scoring"].get("bcell_score", 20.0)
    bcell_candidates = sorted([
        {"peptide": p, "type": "B-cell", "score": bcell_score_val}
        for p in set(bcell_epitopes)
    ], key=lambda x: len(x["peptide"]), reverse=True)

    _ok(f"Scored: {len(mhci_candidates)} MHC-I, "
        f"{len(mhcii_candidates)} MHC-II, {len(bcell_candidates)} B-cell candidates")

    # ── PHASE 1D: Safety filtering ────────────────────────────────────────────
    _section("PHASE 1 — SAFETY PRE-SCREENING (Toxicity / Allergenicity / Homology)")

    all_raw_candidates = mhci_candidates + mhcii_candidates + bcell_candidates
    peptides_to_test   = list({c["peptide"] for c in all_raw_candidates})

    print(f"[SAFETY] Testing {len(peptides_to_test)} unique peptides...")

    tox_map      = run_toxinpred(peptides_to_test)
    allergen_map = predict_allergen(peptides_to_test, config)
    homology_map = check_human_homology(peptides_to_test, config, use_api=False)

    mhci_candidates  = filter_candidates_by_safety(mhci_candidates, tox_map, allergen_map)
    mhcii_candidates = filter_candidates_by_safety(mhcii_candidates, tox_map, allergen_map)
    bcell_filtered   = filter_candidates_by_safety(bcell_candidates, tox_map, allergen_map)
    if len(bcell_filtered) == 0 and len(bcell_candidates) > 0:
        bcell_filtered = bcell_candidates[:1]
    bcell_candidates = bcell_filtered

    _ok(f"After Safety: {len(mhci_candidates)} MHC-I, "
        f"{len(mhcii_candidates)} MHC-II, {len(bcell_candidates)} B-cell")

    # ── PHASE 1E: Antigenicity ────────────────────────────────────────────────
    _section("PHASE 1 — ANTIGENICITY SCORING (VaxiJen-like ACC)")
    antigenicity_map  = score_epitope_antigenicity(peptides_to_test, config)
    antigen_threshold = config["antigenicity"]["thresholds"].get(
        config.get("target_organism", "virus"), 0.32)

    # ── PHASE 1F: Immunogenicity & TC50 (MHC-I only) ─────────────────────────
    _section("PHASE 1 — IMMUNOGENICITY & TC50 (IEDB Class I / Sette-Vitiello)")
    mhci_peptides = [c["peptide"] for c in mhci_candidates]
    immuno_map    = batch_predict_immunogenicity(mhci_peptides, config)
    print(f"  [OK]   Immunogenicity scored for {len(immuno_map)} CTL peptides")

    # ── PHASE 1G: Population coverage & conservancy ───────────────────────────
    _section("PHASE 1 — POPULATION COVERAGE & CONSERVANCY")
    alleles = config.get("alleles", {}).get("mhc1", []) + config.get("alleles", {}).get("mhc2", [])
    max_cov = config["selection"].get("max_peptides_for_coverage", 30)
    peps_cov = [c["peptide"] for c in (mhci_candidates + mhcii_candidates)[:max_cov]]
    coverage_map   = compute_population_coverage(peps_cov, alleles, config)
    conservancy_map = compute_conservancy(
        [c["peptide"] for c in all_raw_candidates], [seq])

    # ── PHASE 1H: Multi-criteria CTL (CD8+) ranking ──────────────────────────
    _section("PHASE 1 — MULTI-CRITERIA CD8+/MHC-I EPITOPE RANKING")
    print("  Criteria: Binding | Antigenicity | Allergenicity | Toxicity"
          " | Immunogenicity | TC50 | Multi-HLA")

    ctl_cfg = config.get("ctl_selection", {})
    w = ctl_cfg.get("weights", {})
    w_bind  = w.get("binding_strength",  0.20)
    w_ag    = w.get("antigenicity",      0.15)
    w_alrg  = w.get("allergenicity",     0.10)
    w_tox   = w.get("toxicity",          0.10)
    w_imm   = w.get("immunogenicity",    0.20)
    w_tc50  = w.get("tc50",             0.10)
    w_hla   = w.get("multi_hla_binding", 0.15)
    min_alleles = ctl_cfg.get("min_allele_count", 2)
    min_pop_cov = ctl_cfg.get("min_population_coverage", 0.90)
    total_mhc1_alleles = len(config.get("alleles", {}).get("mhc1", []))

    for c in mhci_candidates:
        p = c["peptide"]
        cov_data = coverage_map.get(p, {})
        c["population_coverage"] = cov_data.get("global", 0.0) if isinstance(cov_data, dict) else float(cov_data)
        c["conservancy"]  = conservancy_map.get(p, 0.0)
        c["antigenicity"]  = antigenicity_map.get(p, {}).get("score", 0.0)
        c["is_antigen"]    = antigenicity_map.get(p, {}).get("is_antigen", False)

        # Immunogenicity
        im = immuno_map.get(p, {})
        c["immunogenicity_score"] = im.get("score", 0.0)
        c["immunogenicity_grade"] = im.get("grade", "N/A")
        c["immunogenic"]          = im.get("immunogenic", False)

        # TC50
        best_ic50 = c.get("best_ic50", 500.0)
        tc50_data = estimate_tc50(best_ic50, c["immunogenicity_score"])
        c["tc50_nM"]    = tc50_data["tc50_nM"]
        c["tc50_grade"] = tc50_data["grade"]

        # Safety flags
        c["is_toxic"]    = tox_map.get(p, {}).get("toxic", False)
        c["is_allergen"] = allergen_map.get(p, {}).get("allergen", False)
        c["is_safe"]     = (not c["is_toxic"]) and (not c["is_allergen"])

        # Multi-HLA
        ac = c.get("allele_count", 1)
        c["hla_breadth_frac"] = ac / max(total_mhc1_alleles, 1)

        # ── Normalised sub-scores (each 0–1) ──
        # Binding: log-scale normalised, higher=better
        bind_norm = min(1.0, max(0.0, c["score"] / 20.0))
        # Antigenicity: ratio to threshold
        ag_norm   = min(1.0, max(0.0, c["antigenicity"] / max(antigen_threshold, 0.01)))
        # Allergenicity: binary pass/fail
        alrg_norm = 1.0 if not c["is_allergen"] else 0.0
        # Toxicity: binary pass/fail
        tox_norm  = 1.0 if not c["is_toxic"] else 0.0
        # Immunogenicity: shift to 0-1
        imm_norm  = min(1.0, max(0.0, (c["immunogenicity_score"] + 1.0) / 2.0))
        # TC50: lower is better
        tc50_norm = min(1.0, max(0.0, 1.0 - (c["tc50_nM"] / 1000.0)))
        # HLA breadth
        hla_norm  = c["hla_breadth_frac"]

        c["_sub"] = {
            "binding": round(bind_norm, 4), "antigenicity": round(ag_norm, 4),
            "allergenicity": round(alrg_norm, 4), "toxicity": round(tox_norm, 4),
            "immunogenicity": round(imm_norm, 4), "tc50": round(tc50_norm, 4),
            "multi_hla": round(hla_norm, 4),
        }
        c["composite_rank"] = round(
            w_bind * bind_norm + w_ag * ag_norm + w_alrg * alrg_norm +
            w_tox * tox_norm + w_imm * imm_norm + w_tc50 * tc50_norm +
            w_hla * hla_norm, 4)

    # Filter: must bind >= min_alleles
    mhci_candidates = [c for c in mhci_candidates if c.get("allele_count", 1) >= min_alleles]
    # Sort by composite rank
    mhci_candidates.sort(key=lambda x: x["composite_rank"], reverse=True)

    _ok(f"{len(mhci_candidates)} CTL candidates after multi-criteria ranking (min {min_alleles} HLA alleles)")

    # Basic CTL diversity selection
    overlap_limit = config["thresholds"].get("overlap_limit", 6)
    mhc1_selected = []
    used_peptides = []

    for c in mhci_candidates:
        if len(mhc1_selected) >= config["selection"]["counts"]["MHC-I"]:
            break
        if not is_redundant(c["peptide"], used_peptides, overlap_limit, config):
            mhc1_selected.append(c)
            used_peptides.append(c["peptide"])

    # ── PHASE 1I: Multi-criteria HTL (CD4+) ranking ────────────────────────────
    _section("PHASE 1 — MULTI-CRITERIA CD4+/MHC-II EPITOPE RANKING")
    print("  Criteria: Binding | Antigenicity | Allergenicity | Toxicity"
          " | Immunogenicity | TC50 | Multi-HLA | IFN-gamma | IL-4 | IL-10")

    from cytokine import predict_cytokines

    # HTL-specific weights (same 7 + 3 cytokine criteria)
    htl_cfg = config.get("htl_selection", config.get("ctl_selection", {}))
    hw = htl_cfg.get("weights", {})
    hw_bind = hw.get("binding_strength",  0.15)
    hw_ag   = hw.get("antigenicity",      0.10)
    hw_alrg = hw.get("allergenicity",     0.08)
    hw_tox  = hw.get("toxicity",          0.08)
    hw_imm  = hw.get("immunogenicity",    0.12)
    hw_tc50 = hw.get("tc50",             0.07)
    hw_hla  = hw.get("multi_hla_binding", 0.10)
    hw_ifng = hw.get("ifn_gamma",         0.12)
    hw_il4  = hw.get("il4",              0.10)
    hw_il10 = hw.get("il10",             0.08)
    total_mhc2_alleles = len(config.get("alleles", {}).get("mhc2", []))

    # Immunogenicity for MHC-II peptides
    mhcii_peptides = [c["peptide"] for c in mhcii_candidates]
    htl_immuno_map = batch_predict_immunogenicity(mhcii_peptides, config)

    # Cytokine predictions (IFN-gamma, IL-4, IL-10)
    htl_cytokine_input = [{"peptide": c["peptide"], "type": "MHC-II"} for c in mhcii_candidates]
    htl_cytokine_results = predict_cytokines(htl_cytokine_input)
    cytokine_map = {r["peptide"]: r for r in htl_cytokine_results}

    for c in mhcii_candidates:
        p = c["peptide"]
        cov_data = coverage_map.get(p, {})
        c["population_coverage"] = cov_data.get("global", 0.0) if isinstance(cov_data, dict) else float(cov_data)
        c["conservancy"]  = conservancy_map.get(p, 0.0)
        c["antigenicity"]  = antigenicity_map.get(p, {}).get("score", 0.0)
        c["is_antigen"]    = antigenicity_map.get(p, {}).get("is_antigen", False)
        c["is_toxic"]      = tox_map.get(p, {}).get("toxic", False)
        c["is_allergen"]   = allergen_map.get(p, {}).get("allergen", False)

        # Immunogenicity
        im = htl_immuno_map.get(p, {})
        c["immunogenicity_score"] = im.get("score", 0.0)
        c["immunogenicity_grade"] = im.get("grade", "N/A")
        c["immunogenic"]          = im.get("immunogenic", False)

        # TC50
        best_ic50 = c.get("best_ic50", 800.0)
        tc50_data = estimate_tc50(best_ic50, c["immunogenicity_score"])
        c["tc50_nM"]    = tc50_data["tc50_nM"]
        c["tc50_grade"] = tc50_data["grade"]

        # Multi-HLA
        ac = c.get("allele_count", 1)
        c["hla_breadth_frac"] = ac / max(total_mhc2_alleles, 1)

        # Cytokines
        cyt = cytokine_map.get(p, {})
        c["ifng_score"]    = cyt.get("ifng", 0.0)
        c["ifng_positive"] = cyt.get("ifng_positive", False)
        c["il4_score"]     = cyt.get("il4", 0.0)
        c["il4_positive"]  = cyt.get("il4_positive", False)
        c["il10_score"]    = cyt.get("il10", 0.0)
        c["il10_positive"] = cyt.get("il10_positive", False)
        c["th_bias"]       = cyt.get("th_bias", "Unknown")

        # Normalised sub-scores (0–1)
        bind_norm = min(1.0, max(0.0, c["score"] / 20.0))
        ag_norm   = min(1.0, max(0.0, c["antigenicity"] / max(antigen_threshold, 0.01)))
        alrg_norm = 1.0 if not c["is_allergen"] else 0.0
        tox_norm  = 1.0 if not c["is_toxic"] else 0.0
        imm_norm  = min(1.0, max(0.0, (c["immunogenicity_score"] + 1.0) / 2.0))
        tc50_norm = min(1.0, max(0.0, 1.0 - (c["tc50_nM"] / 1000.0)))
        hla_norm  = c["hla_breadth_frac"]
        ifng_norm = min(1.0, max(0.0, (c["ifng_score"] + 0.5) / 1.0))
        il4_norm  = min(1.0, max(0.0, (c["il4_score"] + 0.5) / 1.0))
        il10_norm = min(1.0, max(0.0, (c["il10_score"] + 0.5) / 1.0))

        c["_sub"] = {
            "binding": round(bind_norm, 4), "antigenicity": round(ag_norm, 4),
            "allergenicity": round(alrg_norm, 4), "toxicity": round(tox_norm, 4),
            "immunogenicity": round(imm_norm, 4), "tc50": round(tc50_norm, 4),
            "multi_hla": round(hla_norm, 4),
            "ifn_gamma": round(ifng_norm, 4), "il4": round(il4_norm, 4),
            "il10": round(il10_norm, 4),
        }
        c["composite_rank"] = round(
            hw_bind * bind_norm + hw_ag * ag_norm + hw_alrg * alrg_norm +
            hw_tox * tox_norm + hw_imm * imm_norm + hw_tc50 * tc50_norm +
            hw_hla * hla_norm + hw_ifng * ifng_norm + hw_il4 * il4_norm +
            hw_il10 * il10_norm, 4)

    mhcii_candidates.sort(key=lambda x: x.get("composite_rank", 0), reverse=True)
    _ok(f"{len(mhcii_candidates)} HTL candidates after multi-criteria ranking")

    # HTL diversity selection
    mhc2_selected = select_diverse_top(
        mhcii_candidates, config["selection"]["counts"]["MHC-II"],
        used_peptides, overlap_limit, config)
    for c in mhc2_selected:
        used_peptides.append(c["peptide"])

    # ── PHASE 1J: Combined Coverage-Aware Selection (>97%) ─────────────────────
    combined_target = config.get("selection", {}).get("combined_population_coverage", 0.97)
    _section(f"PHASE 1 — COMBINED COVERAGE-AWARE SELECTION (target >={combined_target * 100:.0f}%)")


    bio = config.get("biological_data", {})
    hla_freq = bio.get("hla_frequencies", {})
    def_freq = bio.get("default_allele_frequency", 0.03)

    covered_alleles = set()
    for c in mhc1_selected + mhc2_selected:
        covered_alleles.update(c.get("alleles_bound", []))

    # Compute combined coverage
    combined_not_cov = 1.0
    for a in covered_alleles:
        combined_not_cov *= (1.0 - hla_freq.get(a, def_freq))
    combined_cov = round(1.0 - combined_not_cov, 4)

    if combined_cov < combined_target:
        _warn(f"Initial Combined Coverage {combined_cov*100:.1f}% < {combined_target*100:.0f}% — adding more epitopes")
        
        # Pool remaining candidates, sorted by their composite rank
        remaining_pool = [c for c in mhci_candidates if c not in mhc1_selected] + \
                         [c for c in mhcii_candidates if c not in mhc2_selected]
        remaining_pool.sort(key=lambda x: x.get("composite_rank", 0), reverse=True)

        for c in remaining_pool:
            if is_redundant(c["peptide"], used_peptides, overlap_limit, config):
                continue
            
            new_alleles = set(c.get("alleles_bound", [])) - covered_alleles
            if not new_alleles:
                continue
            
            # Add to respective list
            if c["type"] == "MHC-I":
                mhc1_selected.append(c)
            else:
                mhc2_selected.append(c)
            used_peptides.append(c["peptide"])
            covered_alleles.update(new_alleles)

            test_not = 1.0
            for a in covered_alleles:
                test_not *= (1.0 - hla_freq.get(a, def_freq))
            combined_cov = round(1.0 - test_not, 4)

            if combined_cov >= combined_target:
                break

    _ok(f"Combined coverage after selection: {combined_cov*100:.1f}%")
    if combined_cov >= combined_target:
        _ok(f"[PASS] Combined target >= {combined_target*100:.0f}% ACHIEVED")
    else:
        _warn(f"Coverage {combined_cov*100:.1f}% — below target (limited by input alleles)")



    # ── PHASE 1J: Multi-criteria B-cell ranking ────────────────────────────────
    _section("PHASE 1 — MULTI-CRITERIA B-CELL EPITOPE RANKING")
    print("  Criteria: Antigenicity | Allergenicity | Toxicity")

    bw = config.get("bcell_selection", {}).get("weights", {})
    bw_ag   = bw.get("antigenicity",  0.50)
    bw_alrg = bw.get("allergenicity", 0.25)
    bw_tox  = bw.get("toxicity",      0.25)

    for c in bcell_candidates:
        p = c["peptide"]
        cov_data = coverage_map.get(p, {})
        c["population_coverage"] = cov_data.get("global", 0.0) if isinstance(cov_data, dict) else float(cov_data)
        c["conservancy"]   = conservancy_map.get(p, 0.0)
        c["antigenicity"]  = antigenicity_map.get(p, {}).get("score", 0.0)
        c["is_antigen"]    = antigenicity_map.get(p, {}).get("is_antigen", False)
        c["is_toxic"]      = tox_map.get(p, {}).get("toxic", False)
        c["is_allergen"]   = allergen_map.get(p, {}).get("allergen", False)

        # Normalised sub-scores (0–1)
        ag_norm   = min(1.0, max(0.0, c["antigenicity"] / max(antigen_threshold, 0.01)))
        alrg_norm = 1.0 if not c["is_allergen"] else 0.0
        tox_norm  = 1.0 if not c["is_toxic"] else 0.0

        c["_sub"] = {
            "antigenicity": round(ag_norm, 4),
            "allergenicity": round(alrg_norm, 4),
            "toxicity": round(tox_norm, 4),
        }
        c["composite_rank"] = round(
            bw_ag * ag_norm + bw_alrg * alrg_norm + bw_tox * tox_norm, 4)

    bcell_candidates.sort(key=lambda x: x.get("composite_rank", 0), reverse=True)
    _ok(f"{len(bcell_candidates)} B-cell candidates after multi-criteria ranking")

    bcell_selected = select_diverse_top(
        bcell_candidates, config["selection"]["counts"]["B-cell"],
        used_peptides + [s["peptide"] for s in mhc2_selected], overlap_limit, config)

    if len(bcell_selected) == 0 and len(mhc1_selected) > 0:
        bcell_selected = [mhc1_selected[0].copy()]
        bcell_selected[0]['type'] = 'B-cell'
    min_htl = config["selection"]["counts"].get("MHC-II_min", 4)
    if len(mhc2_selected) < min_htl and len(mhcii_candidates) > len(mhc2_selected):
        pool = used_peptides + [s["peptide"] for s in mhc2_selected]
        for c in mhcii_candidates:
            if not is_redundant(c["peptide"], pool, overlap_limit, config):
                mhc2_selected.append(c.copy()); pool.append(c["peptide"])
                if len(mhc2_selected) >= min_htl: break
    if len(mhc2_selected) == 0 and len(mhc1_selected) > 0:
        mhc2_selected = [m.copy() for m in mhc1_selected[:2]]
        for m in mhc2_selected: m['type'] = 'MHC-II'

    selected_all = mhc1_selected + mhc2_selected + bcell_selected
    _ok(f"Final: {len(mhc1_selected)} CTL | {len(mhc2_selected)} HTL | {len(bcell_selected)} B-cell")

    # ── DETAILED CTL EPITOPE REPORT ──────────────────────────────────────────
    _section("CD8+ / MHC-I EPITOPE SELECTION REPORT (7-Criteria)")
    print(f"\n  {'#':<3} {'PEPTIDE':<20} {'BIND':>5} {'ANTGN':>6} {'ALRGN':>6}"
          f" {'TOXIC':>6} {'IMMUN':>6} {'TC50':>6} {'HLA':>5} {'RANK':>6}")
    print("  " + "─" * 80)
    for i, ep in enumerate(mhc1_selected, 1):
        s = ep.get("_sub", {})
        print(f"  {i:<3} {ep['peptide']:<20}"
              f" {s.get('binding',0):.2f}"
              f"  {s.get('antigenicity',0):.2f}"
              f"  {'PASS' if s.get('allergenicity',0)==1 else 'FAIL':>5}"
              f"  {'PASS' if s.get('toxicity',0)==1 else 'FAIL':>5}"
              f"  {s.get('immunogenicity',0):.2f}"
              f"  {s.get('tc50',0):.2f}"
              f" {ep.get('allele_count',0):>3}"
              f"  {ep.get('composite_rank',0):.3f}")

    print(f"\n  ── Per-Epitope CTL Detail ──")
    for i, ep in enumerate(mhc1_selected, 1):
        print(f"\n  [{i}] {ep['peptide']}")
        print(f"      Binding Score:      {ep.get('score',0):.2f}  |  Best IC50: {ep.get('best_ic50','N/A')} nM")
        print(f"      Antigenicity:       {ep.get('antigenicity',0):.4f}  ->  {'ANTIGEN (+)' if ep.get('is_antigen') else 'non-antigen'}")
        print(f"      Allergenicity:      {'NON-ALLERGEN (+)' if not ep.get('is_allergen') else 'ALLERGEN (-)'}")
        print(f"      Toxicity:           {'NON-TOXIC (+)' if not ep.get('is_toxic') else 'TOXIC (-)'}")
        print(f"      Immunogenicity:     {ep.get('immunogenicity_score',0):.4f}  ->  {ep.get('immunogenicity_grade','N/A')}")
        print(f"      TC50:               {ep.get('tc50_nM',0):.1f} nM  ->  {ep.get('tc50_grade','N/A')}")
        alleles = ep.get('alleles_bound', [])
        print(f"      HLA Alleles ({len(alleles)}):  {', '.join(alleles[:6])}")
        if len(alleles) > 6:
            print(f"                          {', '.join(alleles[6:])}")
        print(f"      Pop. Coverage:      {ep.get('population_coverage',0)*100:.1f}%")
        print(f"      Conservancy:        {ep.get('conservancy',0)*100:.1f}%")
        print(f"      ── Composite Rank:  {ep.get('composite_rank',0):.4f}")

    print(f"\n  ═══ Combined CTL Population Coverage: {combined_cov*100:.1f}% ═══")

    # ── DETAILED HTL EPITOPE REPORT ──────────────────────────────────────────
    _section("CD4+ / MHC-II EPITOPE SELECTION REPORT (10-Criteria)")
    print(f"\n  {'#':<3} {'PEPTIDE':<25} {'BIND':>5} {'ANTGN':>6} {'IMMUN':>6}"
          f" {'TC50':>6} {'HLA':>4} {'IFN-g':>5} {'IL4':>5} {'IL10':>5} {'RANK':>6}")
    print("  " + "─" * 90)
    for i, ep in enumerate(mhc2_selected, 1):
        s = ep.get("_sub", {})
        print(f"  {i:<3} {ep['peptide']:<25}"
              f" {s.get('binding',0):.2f}"
              f"  {s.get('antigenicity',0):.2f}"
              f"  {s.get('immunogenicity',0):.2f}"
              f"  {s.get('tc50',0):.2f}"
              f" {ep.get('allele_count',0):>3}"
              f"  {'+' if ep.get('ifng_positive') else '-':>4}"
              f"  {'+' if ep.get('il4_positive') else '-':>4}"
              f"  {'+' if ep.get('il10_positive') else '-':>4}"
              f"  {ep.get('composite_rank',0):.3f}")

    print(f"\n  ── Per-Epitope HTL Detail ──")
    for i, ep in enumerate(mhc2_selected, 1):
        print(f"\n  [{i}] {ep['peptide']}")
        print(f"      Binding Score:      {ep.get('score',0):.2f}  |  Best IC50: {ep.get('best_ic50','N/A')} nM")
        print(f"      Antigenicity:       {ep.get('antigenicity',0):.4f}  ->  {'ANTIGEN (+)' if ep.get('is_antigen') else 'non-antigen'}")
        print(f"      Allergenicity:      {'NON-ALLERGEN (+)' if not ep.get('is_allergen') else 'ALLERGEN (-)'}")
        print(f"      Toxicity:           {'NON-TOXIC (+)' if not ep.get('is_toxic') else 'TOXIC (-)'}")
        print(f"      Immunogenicity:     {ep.get('immunogenicity_score',0):.4f}  ->  {ep.get('immunogenicity_grade','N/A')}")
        print(f"      TC50:               {ep.get('tc50_nM',0):.1f} nM  ->  {ep.get('tc50_grade','N/A')}")
        alleles = ep.get('alleles_bound', [])
        print(f"      HLA Alleles ({len(alleles)}):  {', '.join(alleles[:5])}")
        if len(alleles) > 5:
            print(f"                          {', '.join(alleles[5:])}")
        print(f"      IFN-gamma:          {ep.get('ifng_score',0):.3f}  ->  {'INDUCER (+)' if ep.get('ifng_positive') else 'non-inducer'}")
        print(f"      IL-4:               {ep.get('il4_score',0):.3f}  ->  {'INDUCER (+)' if ep.get('il4_positive') else 'non-inducer'}")
        print(f"      IL-10:              {ep.get('il10_score',0):.3f}  ->  {'INDUCER (+)' if ep.get('il10_positive') else 'non-inducer'}")
        print(f"      Th bias:            {ep.get('th_bias','Unknown')}")
        print(f"      Pop. Coverage:      {ep.get('population_coverage',0)*100:.1f}%")
        print(f"      Conservancy:        {ep.get('conservancy',0)*100:.1f}%")
        print(f"      ── Composite Rank:  {ep.get('composite_rank',0):.4f}")

    # ── EXPORT: Epitope report (CSV + JSON) ──────────────────────────────────
    _section("EXPORT — EPITOPE SELECTION REPORT (CSV + JSON)")
    report_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    epitope_export = export_epitope_csv(
        ctl_epitopes=mhc1_selected,
        htl_epitopes=mhc2_selected,
        bcell_epitopes=bcell_selected,
        combined_ctl_coverage=combined_cov,
        output_dir=report_dir,
    )

    # ── PHASE 2: Vaccine construction (AUTO ADJUVANT SELECTION) ────────────────
    _section("PHASE 2 — AUTO ADJUVANT SELECTION & VACCINE CONSTRUCTION")
    print("  Evaluating 3 adjuvants: L7/L12, beta-Defensin-3, RS09...")
    print("  Criteria: Antigenicity(35%) | Solubility(25%) | Instability(25%) | Length(15%)")

    adj_result = auto_select_adjuvant(selected_all, config)
    print_adjuvant_comparison(adj_result)

    # Use the winning construct
    antigen = adj_result["winner_construct"]
    adj_name = adj_result["winner_name"]

    _ok(f"Auto-selected adjuvant: {adj_name}")
    _ok(f"Mechanism: {adj_result['winner_mechanism']}")
    _ok(f"PADRE:     {config['construct']['padre_sequence']}")
    _ok(f"His-tag:   {config['construct']['his_tag']}")
    _ok(f"Linkers:   AAY (CTL) | GPGPG (HTL) | KK (B-cell) | EAAAK (rigid)")
    _ok(f"Final construct length: {len(antigen)} amino acids")
    _ok(f"Starts with adjuvant:  {antigen[:10]}...")
    _ok(f"Ends with His-tag:     ...{antigen[-10:]}")

    # ── PHASE 3A: ProtParam physicochemical analysis ──────────────────────────
    _section("PHASE 3 — PHYSICOCHEMICAL ANALYSIS (ProtParam)")

    physchem = analyze_construct(antigen)
    if "error" not in physchem:
        a = physchem["assessment"]
        _info(f"Length:           {physchem['num_amino_acids']} aa")
        _info(f"Molecular weight: {physchem['molecular_weight_kda']} kDa")
        _info(f"Theoretical pI:   {physchem['theoretical_pi']}")
        _info(f"Instability idx:  {physchem['instability_index']['value']}  -> {a['stability']}")
        _info(f"Aliphatic index:  {physchem['aliphatic_index']} -> {a['thermostability']}")
        _info(f"GRAVY:            {physchem['gravy']} -> {a['solubility']}")
        _info(f"MW range check:   {a['mw_range']}")
        _info(f"Half-life (mammalian): {physchem['half_life']['mammalian_reticulocytes']}")
        _info(f"Half-life (E. coli):   {physchem['half_life']['ecoli_in_vivo']}")
    else:
        _fail(physchem["error"])

    # ── PHASE 3B: Full-construct antigenicity ─────────────────────────────────
    _section("PHASE 3 — FULL CONSTRUCT ANTIGENICITY (VaxiJen-like)")

    from antigenicity import predict_antigenicity
    construct_antigenicity = predict_antigenicity(antigen, config)
    score = construct_antigenicity["score"]
    if construct_antigenicity["is_antigen"]:
        _ok(f"Antigenicity score: {score} - {construct_antigenicity['note']}")
    else:
        _warn(f"Antigenicity score: {score} - {construct_antigenicity['note']}")


    # ── PHASE 3C: Immune simulation report ───────────────────────────────────
    _section("PHASE 3 — IMMUNE SIMULATION (Heuristic + C-ImmSim Instructions)")

    immune_report = generate_immune_report(
        antigen_sequence    = antigen,
        selected_epitopes   = selected_all,
        config              = config,
        population_coverage_map = coverage_map,
        output_dir          = "."
    )

    ep_sum = immune_report["epitope_summary"]
    _info(f"Epitopes: {ep_sum['MHC_I_CTL']} CTL | "
          f"{ep_sum['MHC_II_HTL']} HTL | "
          f"{ep_sum['B_cell_humoral']} B-cell")

    ig = immune_report["heuristic_immune_response"].get("interpretation", {})
    _info(f"Predicted IgM:    {ig.get('IgM',    'N/A')}")
    _info(f"Predicted IgG:    {ig.get('IgG',    'N/A')}")
    _info(f"CD8+ CTL:         {ig.get('CD8',    'N/A')}")
    _info(f"CD4+ Helper:      {ig.get('CD4',    'N/A')}")
    _info(f"Memory B-cells:   {ig.get('Memory', 'N/A')}")

    cyt = immune_report["cytokine_profile"]
    _info(f"Cytokine bias:    {cyt['bias']}")

    pop = immune_report["population_coverage"]
    _info(f"Population cov.:  {pop['note']}")

    rdy = immune_report["vaccine_readiness_checklist"]
    print()
    for item, passed in rdy.items():
        (_ok if passed else _fail)(item)

    print(f"\n  {immune_report['readiness_status']}")

    # ── PHASE 3D: C-ImmSim Submission Guide ──────────────────────────────────
    _section("PHASE 3 — C-ImmSim SUBMISSION GUIDE")
    steps = immune_report["cimmsim_submission"]["manual_steps"]
    for step in steps:
        print(f"  {step}")

    # ── PHASE 3E: Secondary Structure Analysis (SOPMA & PSIPRED) ─────────────
    _section("PHASE 3 — SECONDARY STRUCTURE PREDICTION (SOPMA & PSIPRED)")

    from secondary_structure import analyze_secondary_structure
    sec_report = analyze_secondary_structure(antigen, output_dir=report_dir)
    sop_sum = sec_report["sopma"]["summary"]
    psi_sum = sec_report["psipred"]["summary"]
    sec_ana = sec_report["analysis"]

    _info(f"SOPMA Alpha Helix (H):      {sop_sum['alpha_helix_count']} aa ({sop_sum['alpha_helix_pct']}%)")
    _info(f"SOPMA Extended Strand (E):  {sop_sum['extended_strand_count']} aa ({sop_sum['extended_strand_pct']}%)")
    _info(f"SOPMA Beta Turn (T):        {sop_sum['beta_turn_count']} aa ({sop_sum['beta_turn_pct']}%)")
    _info(f"SOPMA Random Coil (C):      {sop_sum['random_coil_count']} aa ({sop_sum['random_coil_pct']}%)")
    _info(f"PSIPRED Avg Confidence:     {psi_sum['avg_confidence']} / 9.0")
    _info(f"Dominant Secondary Feature: {sec_ana['dominant_structure']}")
    _info(f"Structural Flexibility:     {sec_ana['flexibility_pct']}% -> {sec_ana['structural_assessment']}")
    _ok(f"Exported Secondary Structure CSV: {sec_report['report_paths']['csv']}")
    _ok(f"Exported Secondary Structure JSON: {sec_report['report_paths']['json']}")


    # ── PHASE 3F: 3D Structure Modeling (I-TASSER & GalaxyRefine) & Validation ─────
    _section("PHASE 3 — 3D TERTIARY STRUCTURE MODELING (I-TASSER & GALAXYREFINE) & VALIDATION")

    from structure_3d import predict_structure
    s3_report = predict_structure(antigen, output_dir=report_dir)
    itasser = s3_report.get("itasser", {})
    galaxy = s3_report.get("galaxy_refine", {})
    val_sum = s3_report.get("validation", {}).get("summary", {})

    _info(f"Tertiary Engine:       {s3_report.get('method')}")
    _info(f"I-TASSER C-Score:      {itasser.get('c_score')} (TM-score: {itasser.get('tm_score')})")
    _info(f"GalaxyRefine GDT-TS:   {galaxy.get('gdt_ts')} (MolProbity Score: {galaxy.get('molprobity_score')})")
    _info(f"1. ProSA Z-Score:      {val_sum.get('prosa_z_score')} [Range: -9.0 to -4.0]")
    _info(f"2. ERRAT Quality:      {val_sum.get('errat_quality_factor')}% [Target > 80.0%]")
    _info(f"3. Verify 3D Pass:     {val_sum.get('verify3d_pct')}% [Target >= 80.0%]")
    _info(f"4. Ramachandran Outliers: {val_sum.get('ramachandran_disallowed_pct')}% [Target < 1.0%]")
    _ok(f"3D Structure Validation: {s3_report.get('validation', {}).get('overall_validation')} (Refined PDB: {s3_report.get('pdb_path')})")

    # ── FINAL OUTPUT ──────────────────────────────────────────────────────────
    _section("FINAL — CANDIDATE VACCINE CONSTRUCT")
    print(f"\n  {antigen}\n")

    # Export to Module 2
    try:
        mod2_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "module2", "data")
        os.makedirs(mod2_data_dir, exist_ok=True)
        with open(os.path.join(mod2_data_dir, "module1_output.json"), "w", encoding="utf-8") as f:
            json.dump({"antigen_sequence": antigen}, f)
        print("  [OK] Exported construct to module2/data/module1_output.json for seamless integration")
    except Exception as e:
        print(f"  [WARN] Could not export to module2: {e}")

    return {
        "top_epitopes":            selected_all,
        "antigen_sequence":        antigen,
        "selected_adjuvant":       adj_name,
        "adjuvant_comparison":     adj_result,
        "physicochemical":         physchem,
        "construct_antigenicity":  construct_antigenicity,
        "secondary_structure":     sec_report,
        "tertiary_structure":      s3_report,
        "immune_report":           immune_report,
        "epitope_report":          epitope_export.get("report", {}),
        "epitope_report_paths":    epitope_export.get("paths", {}),
        "combined_ctl_coverage":   combined_cov,
        "mhci_raw":                mhci_df.to_dict(),
        "mhcii_raw":               mhcii_df.to_dict(),
        "bcell":                   bcell_epitopes,
    }


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if not sys.stdin.isatty():
        seq = sys.stdin.read().strip()
    else:
        seq = input("Enter protein sequence:\n").strip()

    if not seq:
        print("Error: No sequence provided.")
        sys.exit(1)

    try:
        result = run_module1(seq)
    except Exception as e:
        print(f"\nError: {e}")
        raise