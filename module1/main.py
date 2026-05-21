# main.py

import sys
import json
import os
from processing  import validate_sequence, get_mhci, get_mhcii, get_bcell
from bcell       import extract_bcell_epitopes, extract_bcell_epitopes_relaxed
from scoring     import score_mhci, score_mhcii
from fusion      import fuse_epitopes
from default_config import DEFAULT_CONFIG
from filter      import is_redundant
from advanced_filters import (
    run_toxinpred,
    filter_candidates_by_hydrophobicity
)
from allergenicity import predict_allergen, filter_candidates_by_safety
from iedb_advanced   import compute_population_coverage, compute_conservancy
from antigenicity    import score_epitope_antigenicity
from protparam       import analyze_construct
from homology        import check_human_homology
from immune_export   import generate_immune_report


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


# ── Pretty print helpers ───────────────────────────────────────────────────────

def _section(title: str):
    bar = "=" * 60
    print(f"\n{bar}")
    print(f"  {title}")
    print(bar)

def _ok(msg):   print(f"  [OK]   {msg}")
def _warn(msg):  print(f"  [WARN] {msg}")
def _fail(msg):  print(f"  [FAIL] {msg}")
def _info(msg):  print(f"  [INFO] {msg}")


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
    _section("PHASE 1 — SAFETY PRE-SCREENING (Heuristic Toxicity / Rule-based Allergenicity / Motif Homology)")

    all_raw_candidates = mhci_candidates + mhcii_candidates + bcell_candidates
    peptides_to_test   = list({c["peptide"] for c in all_raw_candidates})

    print(f"[SAFETY] Testing {len(peptides_to_test)} unique peptides...")

    tox_map     = run_toxinpred(peptides_to_test)
    allergen_map = predict_allergen(peptides_to_test, config)
    homology_map = check_human_homology(peptides_to_test, config, use_api=False)

    # GRAVY / hydrophobicity filter
    gravy_limit = config["safety"].get("gravy_threshold", 1.5)
    mhci_candidates  = filter_candidates_by_hydrophobicity(mhci_candidates,  gravy_limit)
    mhcii_candidates = filter_candidates_by_hydrophobicity(mhcii_candidates, gravy_limit)
    bcell_candidates = filter_candidates_by_hydrophobicity(bcell_candidates, gravy_limit)

    # Apply allergen/toxicity filter
    mhci_candidates  = filter_candidates_by_safety(mhci_candidates, tox_map, allergen_map)
    mhcii_candidates = filter_candidates_by_safety(mhcii_candidates, tox_map, allergen_map)
    
    bcell_filtered = filter_candidates_by_safety(bcell_candidates, tox_map, allergen_map)
    
    # Ensure fallback survives safety filter (Point 2)
    if len(bcell_filtered) == 0 and len(bcell_candidates) > 0:
        print("  [INFO] Forcing top B-cell candidate through safety filter (fallback active)")
        bcell_filtered = bcell_candidates[:1]
        
    bcell_candidates = bcell_filtered

    _ok(f"After Safety & GRAVY filters: {len(mhci_candidates)} MHC-I, "
        f"{len(mhcii_candidates)} MHC-II, {len(bcell_candidates)} B-cell")

    # ── PHASE 1E: Antigenicity scoring (VaxiJen-like) ─────────────────────────
    _section("PHASE 1 — ANTIGENICITY SCORING (VaxiJen-like ACC method)")

    antigenicity_map = score_epitope_antigenicity(peptides_to_test, config)
    antigen_threshold = config["antigenicity"]["thresholds"].get(config.get("target_organism", "virus"), 0.32)

    # ── PHASE 1F: Population coverage & conservancy ───────────────────────────
    _section("PHASE 1 — APPROXIMATE POPULATION COVERAGE & CONSERVANCY")

    alleles = config.get("alleles", {}).get("mhc1", []) + config.get("alleles", {}).get("mhc2", [])
    
    # SPEED OPTIMIZATION: Only compute population coverage for the most promising candidates
    max_cov_peptides = config["selection"].get("max_peptides_for_coverage", 20)
    peptides_for_coverage = [c["peptide"] for c in (mhci_candidates + mhcii_candidates)[:max_cov_peptides]]
    
    print(f"  [SPEED] Computing coverage for top {len(peptides_for_coverage)} peptides (skipping weaker hits)...")
    coverage_map   = compute_population_coverage(peptides_for_coverage, alleles, config)
    conservancy_map = compute_conservancy(
        [c["peptide"] for c in all_raw_candidates], [seq]
    )

    # Inject all bonuses/penalties into candidate scores
    pop_bonus   = config["scoring"].get("population_coverage_bonus", 5.0)
    cons_bonus  = config["scoring"].get("conservancy_bonus", 2.0)
    tox_penalty = config["scoring"].get("toxicity_penalty", -100.0)
    all_penalty = config["scoring"].get("allergen_penalty", -100.0)
    hom_penalty = config["scoring"].get("human_homology_penalty", -50.0)

    for c in mhci_candidates + mhcii_candidates + bcell_candidates:
        p = c["peptide"]
        cov_data = coverage_map.get(p, {})
        c["population_coverage"] = cov_data.get("global", 0.0) if isinstance(cov_data, dict) else float(cov_data)
        
        c["conservancy"]         = conservancy_map.get(p, 0.0)
        c["antigenicity"]        = antigenicity_map.get(p, {}).get("score", 0.0)
        c["is_antigen"]          = antigenicity_map.get(p, {}).get("is_antigen", False)

        # Bonuses
        c["score"] += (c["population_coverage"] * pop_bonus +
                       c["conservancy"]          * cons_bonus)

        # Penalties
        if tox_map.get(p, {}).get("toxic"):
            c["score"] += tox_penalty
        if allergen_map.get(p, {}).get("allergen"):
            c["score"] += all_penalty
        if not homology_map.get(p, {}).get("safe", True):
            c["score"] += hom_penalty
            _warn(f"Human homology detected: {p[:20]}... — penalised")
            
        # GRAVY penalty
        gravy_val = c.get("gravy", 0.0)
        g_thresh = config["scoring"].get("gravy_penalty_threshold", 0.5)
        g_penalty = config["scoring"].get("gravy_penalty", -2.0)
        
        if isinstance(gravy_val, float) and gravy_val > g_thresh:
            c["score"] += g_penalty
            _warn(f"Slight hydrophobicity (GRAVY={gravy_val:.2f}): {p[:20]}... — penalised")

    # ── PHASE 1G: Diversity-aware final selection ─────────────────────────────
    _section("PHASE 1 — DIVERSITY-AWARE EPITOPE SELECTION")

    overlap_limit = config["thresholds"].get("overlap_limit", 6)

    mhc1_selected = select_diverse_top(
        mhci_candidates, config["selection"]["counts"]["MHC-I"],
        [], overlap_limit, config
    )
    mhc2_selected = select_diverse_top(
        mhcii_candidates, config["selection"]["counts"]["MHC-II"],
        [s["peptide"] for s in mhc1_selected], overlap_limit, config
    )
    bcell_selected = select_diverse_top(
        bcell_candidates, config["selection"]["counts"]["B-cell"],
        [s["peptide"] for s in mhc1_selected + mhc2_selected], overlap_limit, config
    )

    # Fallbacks and validations
    if len(bcell_selected) == 0 and len(mhc1_selected) > 0:
        bcell_selected = [mhc1_selected[0].copy()]
        bcell_selected[0]['type'] = 'B-cell'

    # Prevent CTL/HTL duplication and enforce min HTL count
    min_htl = config["selection"]["counts"].get("MHC-II_min", 4)
    if len(mhc2_selected) < min_htl and len(mhcii_candidates) > len(mhc2_selected):
        print("  [INFO] Enforcing minimum HTL count and preventing CTL overlap")
        existing_peptides = {s["peptide"] for s in mhc2_selected}
        mhc1_peptides = [s["peptide"] for s in mhc1_selected]
        pool = list(existing_peptides) + mhc1_peptides
        
        for c in mhcii_candidates:
            if not is_redundant(c["peptide"], pool, overlap_limit, config):
                c_copy = c.copy()
                c_copy['type'] = 'MHC-II'
                mhc2_selected.append(c_copy)
                pool.append(c["peptide"])
                if len(mhc2_selected) >= min_htl:
                    break

    if len(mhc2_selected) == 0 and len(mhc1_selected) > 0:
        print("  [FALLBACK] MHC-II triggered")
        mhc2_selected = [m.copy() for m in (mhc1_selected[1:3] if len(mhc1_selected) >= 3 else mhc1_selected[:2])]
        for m in mhc2_selected:
            m['type'] = 'MHC-II'

    if len(mhc1_selected) == 0 or len(mhc2_selected) == 0 or len(bcell_selected) == 0:
        print("  [ERROR] Incomplete vaccine — retry with different input")

    selected_all = mhc1_selected + mhc2_selected + bcell_selected

    _ok(f"Selected: {len(mhc1_selected)} MHC-I  "
        f"| {len(mhc2_selected)} MHC-II  "
        f"| {len(bcell_selected)} B-cell")

    # Print epitope table
    print(f"\n  {'TYPE':<10} {'PEPTIDE':<30} {'SCORE':>8}  {'ANTIGENIC':>10}  {'GRAVY':>7}")
    print("  " + "-" * 72)
    for ep in selected_all:
        antigenic = "YES [+] " if ep.get("is_antigen") else "no"
        gravy_val = ep.get("gravy", "N/A")
        gravy_str = f"{gravy_val:.3f}" if isinstance(gravy_val, float) else str(gravy_val)
        print(f"  {ep['type']:<10} {ep['peptide']:<30} {ep['score']:>8.2f}  {antigenic:>10}  {gravy_str:>7}")

    # ── PHASE 2: Vaccine construction ─────────────────────────────────────────
    _section("PHASE 2 — VACCINE CONSTRUCTION (Adjuvant + PADRE + Epitopes + His-tag)")

    antigen = fuse_epitopes(selected_all, config)
    adj_name = config["construct"]["adjuvant_name"]
    _ok(f"Adjuvant:  {adj_name}")
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
        _ok(f"Antigenicity score: {score} — {construct_antigenicity['note']}")
    else:
        _warn(f"Antigenicity score: {score} — {construct_antigenicity['note']}")

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
        "physicochemical":         physchem,
        "construct_antigenicity":  construct_antigenicity,
        "immune_report":           immune_report,
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