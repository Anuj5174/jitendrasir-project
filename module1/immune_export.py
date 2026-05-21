# immune_export.py
"""
C-ImmSim immune simulation export module.

Prepares the vaccine construct for submission to the C-ImmSim server:
  https://kraken.iac.rm.cnr.it/C-IMMSIM/index.php

Also generates a local immune response heuristic summary showing
expected antibody types and immune cell activation patterns.
"""

import json
import os
from typing import List, Dict

# ── Heuristic immune response estimators ─────────────────────────────────────

def _estimate_immunoglobulin_response(mhci_count: int, mhcii_count: int,
                                       bcell_count: int, config: dict) -> dict:
    """Heuristic estimate of antibody class dominance based on epitope composition."""
    total = mhci_count + mhcii_count + bcell_count
    if total == 0:
        return {}

    mhci_frac  = mhci_count  / total
    mhcii_frac = mhcii_count / total
    bcell_frac = bcell_count / total

    weights = config["immune_simulation"]["igm_weights"]
    igm_score = round(weights["ctl"] * mhci_frac + weights["htl"] * mhcii_frac + weights["bcell"] * bcell_frac, 2)
    
    ig_w = config["immune_simulation"]["igg_weights"]
    igg_score = round(ig_w["ctl"] * mhci_frac + ig_w["htl"] * mhcii_frac + ig_w["bcell"] * bcell_frac, 2)
    
    cw = config["immune_simulation"]["cell_activation_weights"]
    cd8_score   = round(cw["cd8_mhc1"] * mhci_frac  + cw["cd8_mhc2"] * mhcii_frac, 2)
    cd4_score   = round(cw["cd4_mhc2"] * mhcii_frac + cw["cd4_mhc1"] * mhci_frac, 2)
    memory_bcell = round(cw["memory_mhc2"] * mhcii_frac + cw["memory_bcell"] * bcell_frac, 2)

    thresh = config["immune_simulation"]["response_thresholds"]
    h, m = thresh.get("high", 0.6), thresh.get("moderate", 0.3)

    return {
        "IgM_primary_response":    igm_score,
        "IgG_secondary_response":  igg_score,
        "CD8_CTL_activation":      cd8_score,
        "CD4_helper_activation":   cd4_score,
        "memory_B_cell_induction": memory_bcell,
        "interpretation": {
            "IgM":    "HIGH" if igm_score   > h else "MODERATE" if igm_score   > m else "LOW",
            "IgG":    "HIGH" if igg_score   > h else "MODERATE" if igg_score   > m else "LOW",
            "CD8":    "HIGH" if cd8_score   > h else "MODERATE" if cd8_score   > m else "LOW",
            "CD4":    "HIGH" if cd4_score   > h else "MODERATE" if cd4_score   > m else "LOW",
            "Memory": "HIGH" if memory_bcell > h else "MODERATE" if memory_bcell > m else "LOW",
        }
    }


def _estimate_cytokine_profile(mhci_count: int, mhcii_count: int, config: dict) -> dict:
    """Heuristic cytokine profile based on epitope type balance."""
    total = max(1, mhci_count + mhcii_count)
    mhci_frac  = mhci_count  / total
    mhcii_frac = mhcii_count / total

    weights = config["immune_simulation"]["cytokine_weights"]
    th1_score = round(weights["th1_ctl_contribution"] * mhci_frac + weights["th1_htl_contribution"] * mhcii_frac, 2)
    th2_score = round((1.0 - weights["th1_ctl_contribution"]) * mhci_frac + weights["th2_bcell_contribution"] * mhcii_frac, 2)

    if th1_score > th2_score + 0.2:
        bias = "Th1-dominant (cellular immunity preferred)"
    elif th2_score > th1_score + 0.2:
        bias = "Th2-dominant (humoral immunity preferred)"
    else:
        bias = "Mixed Th1/Th2 (balanced cellular + humoral response)"

    return {
        "Th1_score":  th1_score,
        "Th2_score":  th2_score,
        "IFN_gamma":  "HIGH" if th1_score > 0.5 else "MODERATE",
        "IL_4":       "HIGH" if th2_score > 0.5 else "MODERATE",
        "IL_10":      "MODERATE",
        "TNF_alpha":  "HIGH" if th1_score > 0.5 else "LOW",
        "bias":       bias
    }

# ── C-ImmSim submission preparation ──────────────────────────────────────────

def generate_cimmsim_input(antigen_sequence: str, config: dict) -> dict:
    cimmsim = config.get("cimmsim", {})
    days = cimmsim.get("injection_days", [1, 84, 168])
    url = cimmsim.get("url", "https://kraken.iac.rm.cnr.it/C-IMMSIM/index.php")
    return {
        "server_url": url,
        "tool_description": (
            "C-ImmSim uses a position-specific scoring matrix (PSSM) to simulate "
            "adaptive immunity. It models B-cell, T-cell, and antibody dynamics "
            "over time after repeated antigen injections."
        ),
        "submission_parameters": {
            "antigen_sequence": antigen_sequence,
            "injection_schedule": {f"injection_{i+1}": f"Day {d}" for i, d in enumerate(days)},
            "volume_per_injection_ml": cimmsim.get("injection_volume", 0.5),
            "seed":  12345,
            "steps": cimmsim.get("simulation_steps", 1050),
            "note": "All other parameters can be left at default values"
        },
        "expected_outputs": [
            "Immunoglobulin levels (IgM, IgG1, IgG2) over time",
            "B-cell population dynamics (plasma + memory B cells)",
            "Cytotoxic T-cell (CTL / CD8+) activation",
            "Helper T-cell (CD4+) activation",
            "Antigen presentation by APCs",
        ],
        "interpretation_guide": {
            "IgM_spike":      "Expected after primary injection — marks early immune response",
            "IgG_rise":       "Expected after booster injections — marks memory response",
            "memory_B_cells": "High count = long-lasting immunity",
            "low_antigen":    "Antigen clearance after each injection = healthy response",
        },
        "manual_steps": [
            f"1. Open {url}",
            "2. Paste the vaccine antigen_sequence into 'Antigen Sequence' field",
            f"3. Set injection schedule: {', '.join([f'Day {d}' for d in days])}",
            "4. Click 'Run Simulation' — results appear in ~2 minutes",
            "5. Download and attach the IgG/IgM graph to your report",
        ]
    }


# ── Full immune simulation report ─────────────────────────────────────────────

def generate_immune_report(antigen_sequence: str,
                            selected_epitopes: List[dict],
                            config: dict,
                            population_coverage_map: Dict[str, float] = None,
                            output_dir: str = None) -> dict:
    
    mhci_eps  = [e for e in selected_epitopes if e.get("type") == "MHC-I"]
    mhcii_eps = [e for e in selected_epitopes if e.get("type") == "MHC-II"]
    bcell_eps = [e for e in selected_epitopes if e.get("type") == "B-cell"]

    ig_response   = _estimate_immunoglobulin_response(len(mhci_eps), len(mhcii_eps), len(bcell_eps), config)
    cytokine_prof = _estimate_cytokine_profile(len(mhci_eps), len(mhcii_eps), config)
    cimmsim_data  = generate_cimmsim_input(antigen_sequence, config)

    # Population coverage summary
    if population_coverage_map:
        coverages = []
        for v in population_coverage_map.values():
            cov = v.get("global", 0.0) if isinstance(v, dict) else float(v)
            if cov > 0:
                coverages.append(cov)
        avg_cov = round(sum(coverages) / len(coverages) * 100, 2) if coverages else 0.0
    else:
        avg_cov = 0.0

    construct_cfg = config.get("construct", {})
    adj_prefix = construct_cfg.get("adjuvant_prefix", "")
    his_suffix = construct_cfg.get("his_tag_suffix", "")
    min_len = config["validation"].get("min_construct_length", 100)

    checklist = {
        "Has adjuvant at N-terminus":          antigen_sequence.startswith(adj_prefix) if adj_prefix else True,
        "Has His-tag at C-terminus":           antigen_sequence.endswith(his_suffix) if his_suffix else True,
        "Contains MHC-I CTL epitopes":         len(mhci_eps) > 0,
        "Contains MHC-II HTL epitopes":        len(mhcii_eps) > 0,
        "Contains B-cell epitopes":            len(bcell_eps) > 0,
        f"Construct length in range (>{min_len} aa)": len(antigen_sequence) > min_len,
        "Antigen sequence non-empty":          bool(antigen_sequence.strip()),
    }
    all_passed = all(checklist.values())

    report = {
        "epitope_summary": {
            "total_selected":  len(selected_epitopes),
            "MHC_I_CTL":       len(mhci_eps),
            "MHC_II_HTL":      len(mhcii_eps),
            "B_cell_humoral":  len(bcell_eps),
        },
        "population_coverage": {
            "average_pct":   avg_cov,
            "note": (
                f"Approximate population coverage estimation using HLA frequency heuristics: {avg_cov:.1f}% of global population. "
                "Run IEDB Population Coverage tool for precise values"
            )
        },
        "heuristic_immune_response":   ig_response,
        "cytokine_profile":            cytokine_prof,
        "cimmsim_submission":          cimmsim_data,
        "vaccine_readiness_checklist": checklist,
        "readiness_status": "READY FOR STRUCTURAL VALIDATION [OK]" if all_passed
                            else "[WARN] Construct incomplete — input sequence too short or lacks sufficient epitope diversity for MHC-II/B-cell coverage. Try a longer antigen.",
    }

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, "immune_simulation_report.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        report["saved_to"] = out_path

    return report
