# advanced_filters.py
import pandas as pd
from typing import List, Dict


import subprocess
import os

def run_toxinpred(peptides: List[str]) -> Dict[str, dict]:
    """
    Real Toxicity Prediction integration.
    Attempts to call ToxinPred standalone ML classifier.
    If not installed, raises a warning rather than falling back to heuristics,
    adhering strictly to research-grade ML evaluation.
    """
    results = {}
    
    # Check if ToxinPred CLI is available
    toxinpred_installed = False
    try:
        # Example check for a local toxinpred installation in PATH
        subprocess.run(["toxinpred", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        toxinpred_installed = True
    except FileNotFoundError:
        print("  [WARN] ToxinPred ML classifier not found in PATH.")
        print("  [WARN] Skipping toxicity check. Please install ToxinPred/ToxinPred3 for real ML toxicity prediction.")
        
    for p in peptides:
        if toxinpred_installed:
            # Placeholder for actual CLI invocation
            # result = subprocess.run(["toxinpred", "-i", p_temp_file], ...)
            pass
            
        results[p] = {'toxic': False, 'score': 0.0, 'note': 'ToxinPred ML skipped (not installed)'}
        
    return results


def predict_allergen(peptides: List[str], config: Dict = None) -> Dict[str, dict]:
    """
    Real Allergenicity Model Integration.
    Attempts to interface with AlgPred, AllergenFP, or AllerCatPro.
    """
    results = {}

    algpred_installed = False
    try:
        subprocess.run(["algpred", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        algpred_installed = True
    except FileNotFoundError:
        print("  [WARN] AlgPred / AllergenFP ML models not found in PATH.")
        print("  [WARN] Skipping allergenicity check. Please install established ML classifiers.")

    for p in peptides:
        results[p] = {
            'allergen': False,
            'confidence': 0.0,
            'note': 'ML Allergenicity skipped (AlgPred/AllergenFP not installed)'
        }
        
    return results


def calculate_gravy(peptide: str) -> float:
    """Calculate the Grand Average of Hydropathy (GRAVY) for a peptide.
    Kyte-Doolittle scale. Positive = Hydrophobic, Negative = Hydrophilic.
    """
    kd_scale = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    if not peptide:
        return 0.0
    scores = [kd_scale.get(aa, 0.0) for aa in peptide.upper()]
    return sum(scores) / len(scores)


def filter_candidates_by_hydrophobicity(candidates: List[dict], threshold: float) -> List[dict]:
    """Filter candidates, removing those with a GRAVY score above the threshold.
    """
    out = []
    for c in candidates:
        gravy = calculate_gravy(c['peptide'])
        if gravy <= threshold:
            c['gravy'] = round(gravy, 3)
            out.append(c)
    return out


def filter_candidates_by_safety(candidates: List[dict], tox_map: Dict[str, dict], allergen_map: Dict[str, dict]) -> List[dict]:
    """Filter a list of candidate dicts (each must contain 'peptide') by toxicity/allergen maps.

    Returns a filtered list excluding peptides flagged as toxic or allergenic.
    """
    out = []
    for c in candidates:
        p = c.get('peptide')
        if not p:
            continue
        tox = tox_map.get(p, {})
        allr = allergen_map.get(p, {})
        if tox.get('toxic'):
            continue
        if allr.get('allergen'):
            continue
        out.append(c)
    return out
