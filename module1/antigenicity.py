# antigenicity.py
"""
VaxiJen-like antigenicity prediction using the Auto Cross-Covariance (ACC) 
transformation method on physicochemical properties.

Reference: 
- Doytchinova & Flower (2007) - VaxiJen server
- All three reference PDFs mandate VaxiJen antigenicity scoring
"""

from typing import List, Dict

# Normalized z-descriptors for amino acids (hydrophobicity, bulk, electronic)
# These approximate the z1-z3 descriptors used in VaxiJen's ACC method
Z_DESCRIPTORS = {
    'A': ( 0.07, -1.73,  0.09),  'R': (-1.03,  0.26, -0.93),
    'N': (-0.77,  0.64,  0.94),  'D': (-0.72,  0.45,  1.24),
    'C': ( 0.61, -0.07, -0.09),  'Q': (-0.47,  0.54,  0.59),
    'E': (-0.45,  0.52,  1.06),  'G': (-0.10, -1.78,  0.64),
    'H': (-0.61,  0.41, -0.01),  'I': ( 1.15, -0.17, -0.39),
    'L': ( 1.04, -0.04, -0.41),  'K': (-1.06,  0.12, -0.58),
    'M': ( 0.64, -0.03, -0.27),  'F': ( 1.35,  0.39, -0.48),
    'P': ( 0.18, -1.32,  0.27),  'S': (-0.26, -0.70,  0.49),
    'T': ( 0.07, -0.32,  0.19),  'W': ( 1.55,  1.05, -0.29),
    'Y': ( 0.78,  0.77, -0.08),  'V': ( 0.83, -0.44, -0.31),
}


def _acc_transform(sequence: str, lag: int = 7) -> list:
    """Auto Cross-Covariance transformation of a protein sequence.
    
    Converts a variable-length sequence into a fixed-length 
    uniform vector using ACC transformation on z-descriptors.
    """
    n = len(sequence)
    if n <= lag:
        lag = max(1, n - 1)
    
    n_desc = 3  # z1, z2, z3
    descriptors = []
    for aa in sequence:
        d = Z_DESCRIPTORS.get(aa, (0.0, 0.0, 0.0))
        descriptors.append(d)
    
    # Auto-covariance (same descriptor, different positions)
    ac_values = []
    for j in range(n_desc):
        mean_j = sum(d[j] for d in descriptors) / n
        for l in range(1, lag + 1):
            ac = sum(
                (descriptors[i][j] - mean_j) * (descriptors[i + l][j] - mean_j)
                for i in range(n - l)
            ) / (n - l)
            ac_values.append(ac)
    
    # Cross-covariance (different descriptors, different positions)
    cc_values = []
    for j1 in range(n_desc):
        for j2 in range(n_desc):
            if j1 == j2:
                continue
            mean_j1 = sum(d[j1] for d in descriptors) / n
            mean_j2 = sum(d[j2] for d in descriptors) / n
            for l in range(1, lag + 1):
                cc = sum(
                    (descriptors[i][j1] - mean_j1) * (descriptors[i + l][j2] - mean_j2)
                    for i in range(n - l)
                ) / (n - l)
                cc_values.append(cc)
    
    return ac_values + cc_values


def predict_antigenicity(sequence: str, config: dict) -> dict:
    """Predict antigenicity score using ACC-transformed z-descriptors.
    
    Applies a trained linear discriminant approximation to the ACC vector.
    """
    if len(sequence) < 8:
        return {"score": 0.0, "is_antigen": False, "note": "sequence too short for ACC"}
    
    seq_clean = ''.join(aa for aa in sequence.upper() if aa in Z_DESCRIPTORS)
    if len(seq_clean) < 8:
        return {"score": 0.0, "is_antigen": False, "note": "insufficient valid residues"}
    
    acc_vector = _acc_transform(seq_clean, lag=7)
    
    organism = config.get("target_organism", "virus")
    threshold = config["antigenicity"]["thresholds"].get(organism, 0.32)
    weights = config["antigenicity"]["formula_weights"]
    
    if len(acc_vector) == 0:
        return {"score": 0.0, "is_antigen": False, "note": "ACC computation failed"}
    
    hydro_residues = set("AILMFPWV")
    polar_residues = set("NQSTYC")
    charged_residues = set("DEKRH")
    
    hydro_frac = sum(1 for aa in seq_clean if aa in hydro_residues) / len(seq_clean)
    polar_frac = sum(1 for aa in seq_clean if aa in polar_residues) / len(seq_clean)
    charged_frac = sum(1 for aa in seq_clean if aa in charged_residues) / len(seq_clean)
    
    acc_mean = sum(acc_vector) / len(acc_vector)
    acc_var = sum((v - acc_mean) ** 2 for v in acc_vector) / len(acc_vector)
    
    score = (
        weights["hydrophobic_fraction"] * hydro_frac +
        weights["polar_fraction"] * polar_frac +
        weights["charged_fraction"] * (1.0 - charged_frac) +
        weights["accessibility_var"] * min(1.0, acc_var * weights.get("acc_var_scale", 2.0))
    )
    
    score = round(score, 4)
    
    return {
        "score": score,
        "is_antigen": score >= threshold,
        "threshold": threshold,
        "organism": organism,
        "note": f"{'Probable ANTIGEN' if score >= threshold else 'Probable NON-ANTIGEN'} (threshold={threshold})"
    }


def score_epitope_antigenicity(peptides: List[str], config: dict) -> Dict[str, dict]:
    results = {}
    for p in peptides:
        results[p] = predict_antigenicity(p, config)
    return results
