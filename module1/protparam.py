# protparam.py
"""
ProtParam-like physicochemical analysis of vaccine constructs.

Calculates the properties required by all three reference PDFs:
- Molecular weight
- Theoretical isoelectric point (pI)
- Instability index (must be < 40 for a stable protein)
- Aliphatic index (higher = more thermostable)
- GRAVY (Grand Average of Hydropathicity)
- Amino acid composition
- Extinction coefficient
- Estimated half-life

Reference: 
- Gasteiger et al. (2005) - ExPASy ProtParam
- Martinelli (2022), Mortazavi et al. (2024), Shabbir et al. (2025)
"""

from typing import Dict

# Amino acid molecular weights (Da) - average isotopic
AA_MW = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'Q': 146.15, 'E': 147.13, 'G': 75.03, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 128.17, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
}

# pK values for pI calculation
PK_VALUES = {
    'N_term': 9.69,  'C_term': 2.34,
    'D': 3.65, 'E': 4.25, 'C': 8.18,
    'Y': 10.07, 'H': 6.00, 'K': 10.54, 'R': 12.48,
}

# Kyte-Doolittle hydropathy scale
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Instability index dipeptide weights (DIWV)
# Guruprasad et al. (1990) - subset of critical pairs
DIWV = {
    'WW': 1.0, 'WC': 1.0, 'WM': 24.68, 'WH': 24.68,
    'DG': 1.0, 'DP': -6.54, 'DD': 1.0, 'DE': 2.26,
    'EE': 33.60, 'EK': 1.0, 'ED': -6.54, 'ER': 1.0,
    'KK': 1.0, 'KD': 1.0, 'KE': 1.0, 'KR': 33.60,
    'RR': 58.28, 'RK': 1.0, 'RD': 1.0, 'RE': 1.0,
    'LL': 33.60, 'LG': -14.03, 'LA': 1.0, 'LK': -7.49,
    'FF': 1.0, 'FI': 1.0, 'FY': 33.60, 'FW': 1.0,
    'AA': 1.0, 'AG': 1.0, 'AV': 1.0, 'AL': 1.0,
    'II': 1.0, 'IV': -7.49, 'IA': 1.0, 'IL': 1.0,
    'VV': 1.0, 'VA': 1.0, 'VL': 1.0, 'VG': -7.49,
    'GG': 13.34, 'GA': 1.0, 'GV': 1.0, 'GL': 1.0,
    'PP': 20.26, 'PG': 1.0, 'PA': 20.26, 'PV': 20.26,
    'SS': 1.0, 'ST': 1.0, 'SA': 1.0, 'SG': 1.0,
    'TT': 1.0, 'TS': 1.0, 'TA': 1.0, 'TG': -7.49,
    'NN': 1.0, 'ND': 1.0, 'NR': 1.0, 'NK': 1.0,
    'QQ': 20.26, 'QR': 1.0, 'QK': 1.0, 'QN': 1.0,
    'CC': 1.0, 'CW': 24.68, 'CM': 33.60, 'CY': 1.0,
    'HH': 1.0, 'HW': -1.88, 'HD': 1.0, 'HK': 24.68,
    'MM': -1.88, 'MC': 1.0, 'ML': 1.0, 'MF': 1.0,
    'YY': 13.34, 'YF': 1.0, 'YW': -9.37, 'YR': -15.91,
}


def compute_molecular_weight(sequence: str) -> float:
    """Molecular weight in Daltons. Subtract water for each peptide bond."""
    water = 18.02
    mw = sum(AA_MW.get(aa, 0.0) for aa in sequence)
    mw -= water * (len(sequence) - 1)  # peptide bond formation
    return round(mw, 2)


def compute_amino_acid_composition(sequence: str) -> Dict[str, dict]:
    """Count and percentage of each amino acid."""
    total = len(sequence)
    composition = {}
    for aa in sorted(set(sequence)):
        count = sequence.count(aa)
        composition[aa] = {
            "count": count,
            "percentage": round(count / total * 100, 2)
        }
    return composition


def compute_charged_residues(sequence: str) -> dict:
    """Count positively and negatively charged residues."""
    positive = sum(1 for aa in sequence if aa in 'RKH')
    negative = sum(1 for aa in sequence if aa in 'DE')
    return {"positive": positive, "negative": negative}


def compute_theoretical_pi(sequence: str) -> float:
    """Estimate isoelectric point using Henderson-Hasselbalch iterative method."""
    # Count ionizable residues
    counts = {aa: sequence.count(aa) for aa in 'DECHYKR'}
    
    # Binary search for pH where net charge = 0
    low, high = 0.0, 14.0
    for _ in range(200):  # 200 iterations for precision
        mid = (low + high) / 2.0
        
        # Positive charges
        charge = 0.0
        # N-terminus
        charge += 1.0 / (1.0 + 10 ** (mid - PK_VALUES['N_term']))
        # K, R, H
        for aa, pk in [('K', PK_VALUES['K']), ('R', PK_VALUES['R']), ('H', PK_VALUES['H'])]:
            charge += counts.get(aa, 0) / (1.0 + 10 ** (mid - pk))
        
        # Negative charges
        # C-terminus
        charge -= 1.0 / (1.0 + 10 ** (PK_VALUES['C_term'] - mid))
        # D, E, C, Y
        for aa, pk in [('D', PK_VALUES['D']), ('E', PK_VALUES['E']), 
                       ('C', PK_VALUES['C']), ('Y', PK_VALUES['Y'])]:
            charge -= counts.get(aa, 0) / (1.0 + 10 ** (pk - mid))
        
        if charge > 0:
            low = mid
        else:
            high = mid
    
    return round((low + high) / 2.0, 2)


def compute_instability_index(sequence: str) -> dict:
    """Guruprasad's instability index. < 40 = stable protein."""
    if len(sequence) < 2:
        return {"value": 0.0, "classification": "unknown"}
    
    total = 0.0
    for i in range(len(sequence) - 1):
        dipeptide = sequence[i] + sequence[i + 1]
        total += DIWV.get(dipeptide, 1.0)
    
    ii = (10.0 / len(sequence)) * total
    classification = "stable" if ii < 40 else "unstable"
    
    return {"value": round(ii, 2), "classification": classification}


def compute_aliphatic_index(sequence: str) -> float:
    """Aliphatic index = X(Ala) + 2.9*X(Val) + 3.9*(X(Ile)+X(Leu))
    Higher = more thermostable.
    """
    n = len(sequence)
    if n == 0:
        return 0.0
    
    a = sequence.count('A') / n * 100
    v = sequence.count('V') / n * 100
    i = sequence.count('I') / n * 100
    l = sequence.count('L') / n * 100
    
    return round(a + 2.9 * v + 3.9 * (i + l), 2)


def compute_gravy(sequence: str) -> float:
    """Grand Average of Hydropathicity (GRAVY) using Kyte-Doolittle."""
    if not sequence:
        return 0.0
    scores = [KD_SCALE.get(aa, 0.0) for aa in sequence]
    return round(sum(scores) / len(scores), 3)


def compute_extinction_coefficient(sequence: str) -> dict:
    """Extinction coefficient at 280 nm (M^-1 cm^-1).
    ε = nTyr * 1490 + nTrp * 5500 + nCys * 125 (assuming all Cys form cystines)
    """
    n_tyr = sequence.count('Y')
    n_trp = sequence.count('W')
    n_cys = sequence.count('C')
    
    # All Cys pairs (disulfide)
    ec_cystine = n_tyr * 1490 + n_trp * 5500 + (n_cys // 2) * 125
    # No cystines
    ec_reduced = n_tyr * 1490 + n_trp * 5500
    
    return {
        "all_cystine": ec_cystine,
        "reduced": ec_reduced
    }


def compute_half_life(sequence: str) -> dict:
    """Estimated half-life based on N-end rule (first amino acid)."""
    # N-end rule half-lives
    half_lives = {
        'mammalian': {
            'M': '>30 hours', 'S': '>30 hours', 'A': '>30 hours', 'G': '>30 hours',
            'T': '>30 hours', 'V': '>30 hours', 'P': '>20 hours',
            'I': '~20 hours', 'L': '~5.5 hours', 'F': '~1 hour',
            'Y': '~10 min', 'W': '~3 min', 'K': '~1.3 hours',
            'R': '~1 hour', 'H': '~3.5 hours', 'D': '~1.1 hours',
            'E': '~1 hour', 'N': '~1.4 hours', 'Q': '~0.8 hours', 'C': '~1.2 hours'
        },
        'yeast': {
            'M': '>20 hours', 'S': '>20 hours', 'A': '>20 hours', 'G': '>20 hours',
            'T': '>20 hours', 'V': '>20 hours', 'P': '>20 hours',
            'I': '~30 min', 'L': '~3 min', 'F': '~3 min',
            'Y': '~10 min', 'W': '~3 min', 'K': '~3 min',
            'R': '~2 min', 'H': '~10 min', 'D': '~3 min',
            'E': '~30 min', 'N': '~3 min', 'Q': '~10 min', 'C': '>20 hours'
        },
        'ecoli': {
            'M': '>10 hours', 'S': '>10 hours', 'A': '>10 hours', 'G': '>10 hours',
            'T': '>10 hours', 'V': '>10 hours', 'P': 'N/A',
            'I': '>10 hours', 'L': '~2 min', 'F': '~2 min',
            'Y': '~2 min', 'W': '~2 min', 'K': '~2 min',
            'R': '~2 min', 'H': '>10 hours', 'D': '>10 hours',
            'E': '>10 hours', 'N': '>10 hours', 'Q': '>10 hours', 'C': '>10 hours'
        }
    }
    
    first_aa = sequence[0] if sequence else 'M'
    return {
        "mammalian_reticulocytes": half_lives['mammalian'].get(first_aa, 'N/A'),
        "yeast_in_vivo": half_lives['yeast'].get(first_aa, 'N/A'),
        "ecoli_in_vivo": half_lives['ecoli'].get(first_aa, 'N/A')
    }


def analyze_construct(sequence: str) -> dict:
    """Full ProtParam-like analysis of a vaccine construct.
    
    This is the main function to call after vaccine construction.
    Returns all physicochemical properties required by the PDFs.
    """
    seq = ''.join(aa for aa in sequence.upper() if aa in AA_MW)
    
    if not seq:
        return {"error": "No valid amino acids in sequence"}
    
    mw = compute_molecular_weight(seq)
    pi = compute_theoretical_pi(seq)
    ii = compute_instability_index(seq)
    ai = compute_aliphatic_index(seq)
    gravy = compute_gravy(seq)
    ec = compute_extinction_coefficient(seq)
    hl = compute_half_life(seq)
    charges = compute_charged_residues(seq)
    composition = compute_amino_acid_composition(seq)
    
    return {
        "num_amino_acids": len(seq),
        "molecular_weight_da": mw,
        "molecular_weight_kda": round(mw / 1000, 2),
        "theoretical_pi": pi,
        "instability_index": ii,
        "aliphatic_index": ai,
        "gravy": gravy,
        "extinction_coefficient": ec,
        "half_life": hl,
        "charged_residues": charges,
        "formula_notes": f"{len(seq)} residues, MW={round(mw / 1000, 1)} kDa",
        "assessment": {
            "stability": "STABLE [+]" if ii["classification"] == "stable" else "UNSTABLE [-]",
            "solubility": "HYDROPHILIC (good solubility) [+]" if gravy < 0 else "HYDROPHOBIC (poor solubility) [!]",
            "mw_range": "IN RANGE [+]" if 30000 <= mw <= 110000 else "OUT OF TYPICAL RANGE [!]",
            "thermostability": "HIGH [+]" if ai > 70 else "MODERATE [!]" if ai > 40 else "LOW [-]"
        }
    }
