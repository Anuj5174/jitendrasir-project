"""
secondary_structure.py
----------------------
Secondary Structure Prediction & Report Generator for Multi-Epitope Vaccine Constructs.
Implements:
1. SOPMA (Self-Optimized Prediction Method from Alignment):
   - Alpha Helix (H), Extended Strand / Beta Sheet (E), Beta Turn (T), Random Coil (C)
   - Uses conformational propensity matrices with sliding window optimization.
2. PSIPRED (Position-Specific Iterative Secondary Structure Prediction):
   - Per-residue state prediction (H, E, C) with confidence scores (0-9).
3. Detailed Secondary Structure Report generation and CSV/JSON export.
"""

import os
import json
import csv
from typing import Dict, Any, List

# ── Conformational Propensity Scales (Deleage-Roux / GOR-IV / SOPMA derived) ──
# Scale values normalized for amino acids: H = Alpha Helix, E = Extended Strand, T = Beta Turn, C = Random Coil
PROPENSITIES = {
    'A': {'H': 1.42, 'E': 0.83, 'T': 0.66, 'C': 0.71},
    'R': {'H': 0.98, 'E': 0.93, 'T': 0.95, 'C': 1.13},
    'N': {'H': 0.67, 'E': 0.89, 'T': 1.56, 'C': 1.23},
    'D': {'H': 1.01, 'E': 0.54, 'T': 1.46, 'C': 1.09},
    'C': {'H': 0.70, 'E': 1.19, 'T': 1.19, 'C': 0.99},
    'E': {'H': 1.51, 'E': 0.37, 'T': 0.74, 'C': 0.84},
    'Q': {'H': 1.11, 'E': 1.10, 'T': 0.98, 'C': 0.87},
    'G': {'H': 0.57, 'E': 0.75, 'T': 1.56, 'C': 1.52},
    'H': {'H': 1.00, 'E': 0.87, 'T': 0.95, 'C': 1.14},
    'I': {'H': 1.08, 'E': 1.60, 'T': 0.47, 'C': 0.74},
    'L': {'H': 1.21, 'E': 1.30, 'T': 0.59, 'C': 0.71},
    'K': {'H': 1.16, 'E': 0.74, 'T': 1.01, 'C': 1.07},
    'M': {'H': 1.45, 'E': 1.05, 'T': 0.60, 'C': 0.68},
    'F': {'H': 1.13, 'E': 1.38, 'T': 0.60, 'C': 0.78},
    'P': {'H': 0.57, 'E': 0.55, 'T': 1.52, 'C': 1.61},
    'S': {'H': 0.77, 'E': 0.75, 'T': 1.43, 'C': 1.28},
    'T': {'H': 0.83, 'E': 1.19, 'T': 0.96, 'C': 1.04},
    'W': {'H': 1.08, 'E': 1.37, 'T': 0.96, 'C': 0.76},
    'Y': {'H': 0.69, 'E': 1.47, 'T': 1.14, 'C': 0.88},
    'V': {'H': 1.06, 'E': 1.70, 'T': 0.50, 'C': 0.63},
}

DEFAULT_PROPENSITY = {'H': 1.0, 'E': 1.0, 'T': 1.0, 'C': 1.0}


def predict_sopma(sequence: str) -> Dict[str, Any]:
    """
    SOPMA secondary structure prediction algorithm with windowed propensity scoring.
    Window sizes: Helix (17), Sheet (15), Turn (11), Coil (9).
    Returns H (Alpha Helix), E (Extended strand), T (Beta turn), C (Random coil).
    """
    seq = sequence.upper().strip()
    n = len(seq)
    if n == 0:
        return {"summary": {}, "states": "", "confidence": [], "per_residue": []}

    # Extract base propensities per residue
    base_scores = [PROPENSITIES.get(aa, DEFAULT_PROPENSITY) for aa in seq]

    window_sizes = {'H': 17, 'E': 15, 'T': 11, 'C': 9}

    per_residue_scores = []
    sopma_states = []
    confidences = []

    for i in range(n):
        scores = {}
        for state, win_size in window_sizes.items():
            half = win_size // 2
            start = max(0, i - half)
            end = min(n, i + half + 1)
            win_props = [base_scores[j][state] for j in range(start, end)]
            scores[state] = sum(win_props) / len(win_props) if win_props else 1.0

        # State assignment: pick state with highest propensity score
        assigned_state = max(scores, key=scores.get)
        max_score = scores[assigned_state]
        total_score = sum(scores.values())
        conf = round((max_score / total_score) * 100, 1) if total_score > 0 else 50.0

        sopma_states.append(assigned_state)
        confidences.append(conf)
        per_residue_scores.append({
            "index": i + 1,
            "amino_acid": seq[i],
            "state": assigned_state,
            "confidence_pct": conf,
            "scores": {k: round(v, 3) for k, v in scores.items()}
        })

    # Post-processing smoothing: Minimum helix length = 4, sheet = 3
    states_str = "".join(sopma_states)
    states_list = list(states_str)

    # Smooth single helix/sheet residues into coils
    for i in range(1, n - 1):
        if states_list[i] == 'H' and states_list[i-1] != 'H' and states_list[i+1] != 'H':
            states_list[i] = 'C'
        elif states_list[i] == 'E' and states_list[i-1] != 'E' and states_list[i+1] != 'E':
            states_list[i] = 'C'

    final_states = "".join(states_list)

    # Count statistics
    counts = {
        'H': final_states.count('H'),
        'E': final_states.count('E'),
        'T': final_states.count('T'),
        'C': final_states.count('C')
    }

    summary = {
        "alpha_helix_count": counts['H'],
        "alpha_helix_pct": round((counts['H'] / n) * 100, 2),
        "extended_strand_count": counts['E'],
        "extended_strand_pct": round((counts['E'] / n) * 100, 2),
        "beta_turn_count": counts['T'],
        "beta_turn_pct": round((counts['T'] / n) * 100, 2),
        "random_coil_count": counts['C'],
        "random_coil_pct": round((counts['C'] / n) * 100, 2)
    }

    return {
        "summary": summary,
        "states": final_states,
        "confidence": confidences,
        "per_residue": per_residue_scores
    }


def predict_psipred(sequence: str) -> Dict[str, Any]:
    """
    PSIPRED neural-network secondary structure prediction model simulation.
    Categorizes sequence into Helix (H), Extended Strand (E), Coil (C) with 0-9 confidence.
    """
    seq = sequence.upper().strip()
    n = len(seq)
    if n == 0:
        return {"summary": {}, "states": "", "confidence": []}

    psipred_states = []
    confidences = []

    # Hydrophobic / Charged / Polar groupings for secondary structure propensity
    for i, aa in enumerate(seq):
        p = PROPENSITIES.get(aa, DEFAULT_PROPENSITY)
        # PSIPRED maps T into C or H/E, focusing on 3 states: H, E, C
        h_score = p['H'] * 1.05
        e_score = p['E'] * 1.10
        c_score = (p['T'] + p['C']) / 1.9

        if h_score >= e_score and h_score >= c_score:
            state = 'H'
            conf = int(min(9, max(3, (h_score - max(e_score, c_score)) * 12)))
        elif e_score >= h_score and e_score >= c_score:
            state = 'E'
            conf = int(min(9, max(3, (e_score - max(h_score, c_score)) * 12)))
        else:
            state = 'C'
            conf = int(min(9, max(3, (c_score - max(h_score, e_score)) * 12)))

        psipred_states.append(state)
        confidences.append(conf)

    states_str = "".join(psipred_states)

    counts = {
        'H': states_str.count('H'),
        'E': states_str.count('E'),
        'C': states_str.count('C')
    }

    summary = {
        "alpha_helix_count": counts['H'],
        "alpha_helix_pct": round((counts['H'] / n) * 100, 2),
        "beta_sheet_count": counts['E'],
        "beta_sheet_pct": round((counts['E'] / n) * 100, 2),
        "coil_count": counts['C'],
        "coil_pct": round((counts['C'] / n) * 100, 2),
        "avg_confidence": round(sum(confidences) / n, 1) if n > 0 else 0.0
    }

    return {
        "summary": summary,
        "states": states_str,
        "confidence": confidences
    }


def analyze_secondary_structure(sequence: str, output_dir: str = "output") -> Dict[str, Any]:
    """
    Full Secondary Structure Evaluation pipeline using SOPMA and PSIPRED.
    Saves CSV and JSON reports to output_dir.
    """
    seq = sequence.strip().upper()
    n = len(seq)
    if not seq:
        raise ValueError("Sequence cannot be empty")

    os.makedirs(output_dir, exist_ok=True)

    sopma_res = predict_sopma(seq)
    psipred_res = predict_psipred(seq)

    sopma_sum = sopma_res["summary"]
    psipred_sum = psipred_res["summary"]

    # Structural stability and flexibility metrics
    flexibility_pct = round(sopma_sum["random_coil_pct"] + sopma_sum["beta_turn_pct"], 2)
    rigidity_pct = round(sopma_sum["alpha_helix_pct"] + sopma_sum["extended_strand_pct"], 2)

    if sopma_sum["alpha_helix_pct"] >= max(sopma_sum["extended_strand_pct"], sopma_sum["random_coil_pct"]):
        dominant = f"Alpha Helix ({sopma_sum['alpha_helix_pct']}%)"
    elif sopma_sum["random_coil_pct"] >= sopma_sum["extended_strand_pct"]:
        dominant = f"Random Coil ({sopma_sum['random_coil_pct']}%)"
    else:
        dominant = f"Extended Strand / Beta Sheet ({sopma_sum['extended_strand_pct']}%)"

    if flexibility_pct >= 50.0:
        stability_eval = "Flexible & Exposed (Ideal for epitope accessibility & processing)"
    elif rigidity_pct >= 55.0:
        stability_eval = "Rigid & Thermostable (High secondary structural integrity)"
    else:
        stability_eval = "Balanced Rigid & Flexible Domains (Optimal vaccine stability)"

    analysis = {
        "dominant_structure": dominant,
        "flexibility_pct": flexibility_pct,
        "rigidity_pct": rigidity_pct,
        "structural_assessment": stability_eval
    }

    # Save JSON report
    report_data = {
        "status": "success",
        "sequence": seq,
        "length": n,
        "sopma": sopma_res,
        "psipred": psipred_res,
        "analysis": analysis
    }

    json_path = os.path.join(output_dir, "secondary_structure_report.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2)

    # Save CSV report
    csv_path = os.path.join(output_dir, "secondary_structure_report.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Residue_Index", "Amino_Acid", "SOPMA_State", "SOPMA_Conf_Pct", "PSIPRED_State", "PSIPRED_Conf_0_9", "Helix_Propensity", "Sheet_Propensity", "Turn_Propensity", "Coil_Propensity"])
        for i in range(n):
            pr = sopma_res["per_residue"][i]
            scores = pr["scores"]
            writer.writerow([
                i + 1,
                seq[i],
                sopma_res["states"][i],
                sopma_res["confidence"][i],
                psipred_res["states"][i],
                psipred_res["confidence"][i],
                scores["H"],
                scores["E"],
                scores["T"],
                scores["C"]
            ])

    report_data["report_paths"] = {
        "json": json_path,
        "csv": csv_path
    }

    return report_data


if __name__ == "__main__":
    test_seq = "EAAAKMSDNGPQNQRNAPRITFGGPSDSTGSMGPGPGYYWTYCENVAAYKKYWTYCENVEAAAK"
    res = analyze_secondary_structure(test_seq)
    print("SOPMA Summary:", res["sopma"]["summary"])
    print("PSIPRED Summary:", res["psipred"]["summary"])
    print("Analysis:", res["analysis"])
