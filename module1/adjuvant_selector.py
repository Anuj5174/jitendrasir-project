"""
adjuvant_selector.py — Automatic adjuvant selection for vaccine constructs.

Given the selected epitopes from the pipeline, evaluates all 3 candidate
adjuvants (L7/L12, beta-Defensin-3, RS09) and returns the best one based on:
  1. Antigenicity      (35%)
  2. Solubility/GRAVY  (25%)
  3. Instability Index (25%)
  4. Construct Length   (15%)

Called automatically by main.py during vaccine construction phase.
"""

from protparam import analyze_construct
from antigenicity import predict_antigenicity

# ── Three candidate adjuvants ─────────────────────────────────────────────
CANDIDATE_ADJUVANTS = {
    "50S Ribosomal Protein L7/L12 (M. tuberculosis)": {
        "sequence": (
            "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVE"
            "AAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKE"
            "AADEAKAKLEAAGATVTVK"
        ),
        "mechanism": "TLR4 agonist, strong Th1 response, proven in TB vaccines",
    },
    "beta-Defensin-3 (hBD-3)": {
        "sequence": "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK",
        "mechanism": "TLR1/2 agonist, chemoattractant for DCs, bridges innate-adaptive",
    },
    "RS09 (Synthetic TLR4 agonist)": {
        "sequence": "APPHALS",
        "mechanism": "Minimal synthetic TLR4 peptide, compact, low immunogenic noise",
    },
}

# Evaluation weights (user-defined priorities)
EVAL_WEIGHTS = {
    "antigenicity":     0.35,
    "solubility":       0.25,
    "instability":      0.25,
    "construct_length": 0.15,
}


def _build_test_construct(adj_seq, epitopes, config):
    """Build a construct using fusion.py with a specific adjuvant sequence."""
    from fusion import fuse_epitopes
    # Temporarily swap adjuvant in config
    original_adj = config["construct"]["adjuvant_sequence"]
    config["construct"]["adjuvant_sequence"] = adj_seq
    construct = fuse_epitopes(epitopes, config)
    config["construct"]["adjuvant_sequence"] = original_adj  # restore
    return construct


def _solubility_label(gravy):
    if gravy < -0.5:  return "HIGH solubility"
    if gravy < 0:     return "GOOD solubility"
    if gravy < 0.5:   return "MODERATE solubility"
    return "POOR solubility"


def auto_select_adjuvant(selected_epitopes, config):
    """
    Evaluate all 3 adjuvants with the given epitopes, rank them,
    and return the best adjuvant's details + construct.

    Parameters:
        selected_epitopes: list of dicts with 'peptide' and 'type' keys
        config: the pipeline config dict

    Returns:
        dict with keys:
          - winner_name, winner_sequence, winner_construct, winner_score
          - all_results: full comparison data
          - ranking: ordered list of (name, score)
    """
    results = {}

    for name, info in CANDIDATE_ADJUVANTS.items():
        adj_seq = info["sequence"]
        construct = _build_test_construct(adj_seq, selected_epitopes, config)

        pp = analyze_construct(construct)
        ag = predict_antigenicity(construct, config)
        ii = pp["instability_index"]

        results[name] = {
            "adjuvant_sequence": adj_seq,
            "mechanism": info["mechanism"],
            "construct": construct,
            "sequence_length": len(construct),
            "molecular_weight_kda": pp["molecular_weight_kda"],
            "theoretical_pi": pp["theoretical_pi"],
            "antigenicity_score": ag["score"],
            "is_antigenic": ag["is_antigen"],
            "instability_index": ii["value"],
            "stability": ii["classification"],
            "gravy": pp["gravy"],
            "solubility": _solubility_label(pp["gravy"]),
            "aliphatic_index": pp["aliphatic_index"],
            "extinction_coefficient": pp["extinction_coefficient"],
            "half_life": pp["half_life"],
            "charged_residues": pp["charged_residues"],
        }

    # ── Normalized weighted scoring ──────────────────────────────────────
    names = list(results.keys())
    scores = {}

    for n in names:
        r = results[n]

        # Antigenicity: higher is better
        ag_vals = [results[x]["antigenicity_score"] for x in names]
        ag_range = max(ag_vals) - min(ag_vals)
        ag_norm = (r["antigenicity_score"] - min(ag_vals)) / (ag_range + 1e-9) if ag_range > 0 else 0.5

        # Solubility: more negative GRAVY = better
        gr_vals = [results[x]["gravy"] for x in names]
        gr_range = max(gr_vals) - min(gr_vals)
        sol_norm = (max(gr_vals) - r["gravy"]) / (gr_range + 1e-9) if gr_range > 0 else 0.5

        # Instability: lower is better
        ii_vals = [results[x]["instability_index"] for x in names]
        ii_range = max(ii_vals) - min(ii_vals)
        ii_norm = (max(ii_vals) - r["instability_index"]) / (ii_range + 1e-9) if ii_range > 0 else 0.5

        # Construct length: 200-400 ideal
        length = r["sequence_length"]
        if 200 <= length <= 400:
            len_norm = 1.0
        elif length < 200:
            len_norm = max(0, length / 200)
        else:
            len_norm = max(0, 1.0 - (length - 400) / 400)

        composite = (
            EVAL_WEIGHTS["antigenicity"]     * ag_norm +
            EVAL_WEIGHTS["solubility"]       * sol_norm +
            EVAL_WEIGHTS["instability"]      * ii_norm +
            EVAL_WEIGHTS["construct_length"] * len_norm
        )
        scores[n] = round(composite, 4)
        results[n]["composite_score"] = scores[n]

    ranking = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    winner_name = ranking[0][0]
    winner = results[winner_name]

    return {
        "winner_name": winner_name,
        "winner_sequence": winner["adjuvant_sequence"],
        "winner_construct": winner["construct"],
        "winner_score": winner["composite_score"],
        "winner_mechanism": winner["mechanism"],
        "all_results": results,
        "ranking": ranking,
        "weights": EVAL_WEIGHTS,
    }


def print_adjuvant_comparison(selection_result):
    """Pretty-print the adjuvant comparison report."""
    ranking = selection_result["ranking"]
    results = selection_result["all_results"]
    weights = selection_result["weights"]

    print(f"\n  Evaluation Weights: Antigenicity={weights['antigenicity']:.0%}  "
          f"Solubility={weights['solubility']:.0%}  "
          f"Instability={weights['instability']:.0%}  "
          f"Length={weights['construct_length']:.0%}\n")

    for i, (name, score) in enumerate(ranking, 1):
        r = results[name]
        tag = " << BEST" if i == 1 else ""
        print(f"  {'-' * 80}")
        print(f"  #{i}  {name}{tag}")
        print(f"  {'-' * 80}")
        print(f"      Composite Score:    {score:.4f}")
        print(f"      Length:             {r['sequence_length']} aa  |  MW: {r['molecular_weight_kda']} kDa")
        print(f"      Antigenicity:       {r['antigenicity_score']:.4f}  ({'ANTIGENIC' if r['is_antigenic'] else 'NON-ANTIGENIC'})")
        print(f"      Instability Index:  {r['instability_index']:.2f}  ({r['stability'].upper()})")
        print(f"      GRAVY:              {r['gravy']:.3f}  ({r['solubility']})")
        print(f"      Aliphatic Index:    {r['aliphatic_index']:.2f}")
        print(f"      Mechanism:          {r['mechanism']}")
        print()

    print(f"  ==> AUTO-SELECTED: {selection_result['winner_name']}")
    print(f"      Score: {selection_result['winner_score']:.4f}")
    print()
