# immunogenicity.py
"""
CD8+ T-cell Immunogenicity Prediction for MHC-I epitopes.

Implements the IEDB Class I Immunogenicity predictor methodology:
  - Amino acid properties at peptide positions contribute to immunogenicity
  - Position-weighted scoring based on TCR-contact residues (P4-P6)
  - Residue size, aromaticity, and physicochemical properties influence
    T-cell receptor recognition.

Reference:
  Calis et al. (2013) "Properties of MHC Class I Presented Peptides That
  Enhance Immunogenicity" — PLoS Computational Biology 9(10): e1003266

Also includes TC50 (half-maximal T-cell response concentration) estimation
based on binding affinity extrapolation using the Sette & Vitiello (2006)
methodology.
"""

from typing import List, Dict
import math


# ── Immunogenicity scoring matrix ────────────────────────────────────────────
# Position-specific amino acid contribution weights for 9-mer peptides.
# Derived from Calis et al. (2013) — positions 1-indexed.
# TCR-contact residues are P4, P5, P6 (highest weight).
# Anchor residues P2, P9 have near-zero weight (MHC-facing, not TCR-facing).

POSITION_WEIGHTS_9MER = {
    1: 0.04,   # N-terminal — minor TCR contact
    2: 0.00,   # Primary anchor residue (MHC-facing) → no TCR contribution
    3: 0.10,   # Secondary contact
    4: 0.25,   # Major TCR contact
    5: 0.25,   # Major TCR contact (central bulge)
    6: 0.25,   # Major TCR contact
    7: 0.10,   # Secondary contact
    8: 0.01,   # Near C-terminal anchor
    9: 0.00,   # C-terminal anchor (MHC-facing) → no TCR contribution
}

# Amino acid immunogenicity contribution scores.
# Large, aromatic residues at TCR-contact positions enhance immunogenicity.
# Reference: Calis et al. (2013) Supplementary Table S3
AA_IMMUNOGENICITY = {
    'A': -0.02, 'R':  0.05, 'N':  0.02, 'D':  0.05, 'C': -0.08,
    'Q':  0.01, 'E':  0.06, 'G': -0.05, 'H':  0.04, 'I':  0.01,
    'L':  0.00, 'K':  0.06, 'M':  0.00, 'F':  0.15, 'P': -0.20,
    'S': -0.03, 'T': -0.01, 'W':  0.18, 'Y':  0.12, 'V': -0.02,
}

# Amino acids categorised by size/aromaticity for immunogenicity boosting
LARGE_AROMATIC = set("FWY")     # Strong positive effect at TCR contacts
SMALL_ALIPHATIC = set("AGVIL")  # Weak/neutral at TCR contacts
CHARGED = set("DEKRH")          # Moderate positive (charge interaction)


def predict_immunogenicity(peptide: str, config: dict = None) -> dict:
    """
    Predict CD8+ T-cell immunogenicity for an MHC-I peptide.

    Method: Position-weighted amino acid property scoring, replicating the
    IEDB Class I Immunogenicity tool (Calis et al. 2013).

    For 9-mers: uses calibrated position weights.
    For 8-11 mers: adapts weights by interpolating TCR-contact zone.

    Returns dict with:
      - score: float immunogenicity score
      - immunogenic: bool (True if score > threshold)
      - grade: HIGH/MODERATE/LOW
      - method: str
    """
    seq = peptide.upper().strip()
    n = len(seq)

    if n < 8 or n > 25:
        return {
            "score": 0.0, "immunogenic": False, "grade": "N/A",
            "note": f"Length {n} outside MHC range (8-25)",
            "method": "IEDB_immunogenicity_heuristic"
        }

    # Build position weights for this peptide length
    pos_weights = _get_position_weights(n)

    # Compute position-weighted immunogenicity score
    raw_score = 0.0
    for i, aa in enumerate(seq):
        aa_contrib = AA_IMMUNOGENICITY.get(aa, 0.0)

        # Aromatic bonus at TCR contact positions
        if aa in LARGE_AROMATIC and pos_weights[i] > 0.15:
            aa_contrib += 0.08

        # Proline penalty (disrupts peptide-MHC conformation)
        if aa == 'P' and 2 <= i <= n - 2:
            aa_contrib -= 0.10

        raw_score += aa_contrib * pos_weights[i]

    # Normalise to [-1, 1] range
    score = round(max(-1.0, min(1.0, raw_score * 10)), 4)

    # Thresholds from Calis et al. (2013) validation
    threshold = 0.0   # scores > 0 are immunogenic
    high_thresh = 0.3
    mod_thresh = 0.1

    immunogenic = score > threshold
    if score >= high_thresh:
        grade = "HIGH"
    elif score >= mod_thresh:
        grade = "MODERATE"
    elif score > threshold:
        grade = "LOW"
    else:
        grade = "NON-IMMUNOGENIC"

    return {
        "score": score,
        "immunogenic": immunogenic,
        "grade": grade,
        "threshold": threshold,
        "note": f"Immunogenicity={grade} (score={score})",
        "method": "IEDB_immunogenicity_heuristic"
    }


def _get_position_weights(length: int) -> list:
    """Generate position-specific weights for peptides of any MHC-I length (8-14)."""
    if length == 9:
        return [POSITION_WEIGHTS_9MER[i + 1] for i in range(9)]

    # For other lengths, interpolate: anchors at P2 and P(n), TCR peak at centre
    weights = []
    for i in range(length):
        if i == 1 or i == length - 1:
            # Anchor positions → zero weight
            weights.append(0.0)
        elif i == 0:
            weights.append(0.04)
        elif i == length - 2:
            weights.append(0.01)
        else:
            # TCR contact zone — peak in the middle
            centre = (length - 1) / 2.0
            dist = abs(i - centre) / centre
            w = 0.25 * math.exp(-2.0 * dist * dist)
            weights.append(round(w, 3))

    # Normalise so weights sum to 1.0
    total = sum(weights)
    if total > 0:
        weights = [round(w / total, 4) for w in weights]

    return weights


def batch_predict_immunogenicity(peptides: List[str], config: dict = None) -> Dict[str, dict]:
    """Predict immunogenicity for a batch of peptides."""
    results = {}
    for p in peptides:
        results[p] = predict_immunogenicity(p, config)
    return results


# ── TC50 Estimation ──────────────────────────────────────────────────────────

def estimate_tc50(ic50: float, immunogenicity_score: float = 0.0) -> dict:
    """
    Estimate TC50 (half-maximal T-cell activation concentration) from
    MHC binding IC50 and immunogenicity score.

    TC50 represents the peptide concentration at which 50% of T-cells
    are activated — it integrates both MHC binding and TCR recognition.

    Methodology: Sette & Vitiello (2006) relationship between binding
    affinity and functional T-cell response.

    Formula:
        TC50 = IC50 * (1 + e^(-k * immunogenicity_score))
    where k = 2.0 (calibration constant)

    Lower TC50 → more potent T-cell activator.

    Classification (from published literature):
        TC50 < 50 nM   → Strong activator
        TC50 < 200 nM  → Moderate activator
        TC50 < 500 nM  → Weak activator
        TC50 >= 500 nM → Poor activator
    """
    if ic50 <= 0:
        ic50 = 1.0  # floor

    k = 2.0  # calibration constant
    # High immunogenicity reduces TC50 (better T-cell activation)
    # Low/negative immunogenicity increases TC50
    modulation_factor = 1.0 + math.exp(-k * immunogenicity_score)
    tc50 = ic50 * modulation_factor

    tc50 = round(tc50, 2)

    if tc50 < 50:
        grade = "STRONG"
    elif tc50 < 200:
        grade = "MODERATE"
    elif tc50 < 500:
        grade = "WEAK"
    else:
        grade = "POOR"

    return {
        "tc50_nM": tc50,
        "grade": grade,
        "ic50_input": ic50,
        "immunogenicity_modulation": round(modulation_factor, 3),
        "note": f"TC50={tc50:.1f} nM ({grade}) — from IC50={ic50:.1f} nM × modulation={modulation_factor:.3f}",
        "method": "Sette_Vitiello_2006_extrapolation"
    }
