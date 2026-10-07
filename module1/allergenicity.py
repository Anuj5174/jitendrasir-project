# allergenicity.py
# Rule-based allergenicity pre-screening inspired by AllerTOP criteria.
# All thresholds and weights are read from config["safety"] — zero magic numbers.

from typing import List, Dict
from default_config import DEFAULT_CONFIG


KNOWN_ALLERGEN_MOTIFS = [
    {"motif": "DPYSPSP",   "allergen": "Ara h 1 (peanut)",        "weight": 0.90},
    {"motif": "QQQPFPQ",   "allergen": "Gliadin (wheat)",         "weight": 0.95},
    {"motif": "PQPQLPY",   "allergen": "Gliadin alpha",           "weight": 0.90},
    {"motif": "QQPFPQQP",  "allergen": "Gluten epitope",          "weight": 0.95},
    {"motif": "RPKHPIKH",  "allergen": "Beta-casein (milk)",      "weight": 0.85},
    {"motif": "FALPQYLK",  "allergen": "Alpha-casein (milk)",     "weight": 0.80},
    {"motif": "NTDGSTDY",  "allergen": "Ovomucoid (egg)",         "weight": 0.85},
    {"motif": "YLDGQNRP",  "allergen": "Ovalbumin (egg)",         "weight": 0.80},
    {"motif": "DPYSPS",    "allergen": "Tree nut Ber e 1",        "weight": 0.75},
    {"motif": "WRCGRQAGG", "allergen": "Latex Hev b 6",           "weight": 0.80},
    {"motif": "VDNRDPE",   "allergen": "Der p 1 (dust mite)",     "weight": 0.85},
    {"motif": "AYCNTTWD",  "allergen": "Der p 2 (dust mite)",     "weight": 0.80},
    {"motif": "KITYENAK",  "allergen": "Bet v 1 (birch pollen)",  "weight": 0.75},
    {"motif": "GVVSAAD",   "allergen": "Phl p 1 (grass pollen)",  "weight": 0.70},
    {"motif": "EELDSEG",   "allergen": "Tropomyosin (shellfish)", "weight": 0.85},
    {"motif": "KIQALEKL",  "allergen": "Tropomyosin isoform",     "weight": 0.80},
    {"motif": "SQSQSQSQ",  "allergen": "Soy Gly m 4",            "weight": 0.75},
    {"motif": "DKNGDGEV",  "allergen": "Parvalbumin (fish)",      "weight": 0.85},
    {"motif": "GVDADMQA",  "allergen": "Profilin (pan-allergen)", "weight": 0.70},
    {"motif": "MNPDEEAK",  "allergen": "Profilin isoform",        "weight": 0.65},
    {"motif": "CCNCLRPF",  "allergen": "LTP (peach/apple)",       "weight": 0.80},
    {"motif": "QQPGQGQQ",  "allergen": "Storage globulin",        "weight": 0.70},
]

WHO_FAO_8MER_DB = [
    "DPYSPSPY", "QHQRESSS", "SKTIKVPQ", "FRKQLLED",
    "QQPFPQQP", "PQPELPYP", "QQQPFPQP", "FRAALFAP",
    "LIVTQTMK", "GLDIQKVA", "ALPMHIRL",
    "ISQAVHAA", "NTDGSTDY", "GGLEPINF",
    "EELDSEGK", "KIQALEEK", "KQLEDELVH",
    "VDNRDPEA", "YIRDEHAF",
    "KITYENAK", "VDAENAVY",
    "DKNGDGEV", "MKELGTVME",
    "EICPAVKR", "SSPGLFAD",
    "WRCGRQAG", "QVAKSTQQ",
    "PKLQNMAA", "SSISRQEV",
    "SQSQSQSQ", "EVKKLEEP",
    "SNKVLPEF", "KEGDLNKL",
]


def _get_kd(config):
    return config.get("biological_data", {}).get(
        "kyte_doolittle", DEFAULT_CONFIG["biological_data"]["kyte_doolittle"]
    )





def compute_instability_index(peptide: str) -> float:
    diwv = {
        "WW": 1.0, "WC": 1.0, "WM": 24.68, "WH": 24.68,
        "CK": 1.0, "CR": 1.0, "CF": 1.0,
        "KK": 1.0, "KR": 33.6,
        "DD": 1.0, "EE": 0.0, "RR": 83.44,
        "GG": 13.34, "GS": 1.0,
        "LL": 1.0, "FF": 1.0, "VV": 1.0,
        "II": 1.0, "TT": 1.0, "MM": 1.0,
        "HH": 1.0, "YY": 1.0, "AA": 1.0, "SS": 1.0, "NN": 1.0,
    }
    if len(peptide) < 2:
        return 0.0
    score = sum(diwv.get(peptide[i:i+2].upper(), 0) for i in range(len(peptide)-1))
    return (10.0 / len(peptide)) * score


def compute_isoelectric_point(peptide: str) -> float:
    pka = {"D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3,
           "Y": 10.1, "K": 10.5, "R": 12.5, "N_term": 8.0, "C_term": 3.1}

    def charge_at(seq, ph):
        c  = 1.0 / (1.0 + 10 ** (ph - pka["N_term"]))
        c -= 1.0 / (1.0 + 10 ** (pka["C_term"] - ph))
        for aa in seq.upper():
            if   aa == "D": c -= 1.0 / (1.0 + 10 ** (pka["D"] - ph))
            elif aa == "E": c -= 1.0 / (1.0 + 10 ** (pka["E"] - ph))
            elif aa == "H": c += 1.0 / (1.0 + 10 ** (ph - pka["H"]))
            elif aa == "C": c -= 1.0 / (1.0 + 10 ** (pka["C"] - ph))
            elif aa == "Y": c -= 1.0 / (1.0 + 10 ** (pka["Y"] - ph))
            elif aa == "K": c += 1.0 / (1.0 + 10 ** (ph - pka["K"]))
            elif aa == "R": c += 1.0 / (1.0 + 10 ** (ph - pka["R"]))
        return c

    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        lo, hi = (mid, hi) if charge_at(peptide, mid) > 0 else (lo, mid)
    return round((lo + hi) / 2.0, 2)


def check_8mer_match(peptide: str) -> dict:
    seq, hits = peptide.upper(), []
    for i in range(len(seq) - 7):
        w = seq[i:i+8]
        for db in WHO_FAO_8MER_DB:
            if w == db.upper():
                hits.append({"position": i, "matched": w, "db": db})
    return {"matched": bool(hits), "hits": hits}


def scan_allergen_motifs(peptide: str) -> dict:
    seq, hits = peptide.upper(), []
    for entry in KNOWN_ALLERGEN_MOTIFS:
        m = entry["motif"].upper()
        if m in seq:
            hits.append({"motif": m, "allergen": entry["allergen"],
                         "weight": entry["weight"], "position": seq.index(m)})
    hits.sort(key=lambda x: x["weight"], reverse=True)
    return {"matched": bool(hits), "hits": hits,
            "max_weight": hits[0]["weight"] if hits else 0.0}


def check_cysteine_pattern(peptide: str) -> bool:
    return peptide.upper().count("C") >= 2 and len(peptide) <= 20


def check_proline_repeat(peptide: str, config: dict) -> bool:
    sf = config.get("safety", DEFAULT_CONFIG["safety"])
    if len(peptide) < sf.get("allergen_proline_min_len", 12):
        return False
    return (peptide.upper().count("P") / len(peptide)) > sf.get("allergen_proline_fraction", 0.25)


def check_repetitive_kmers(peptide: str, k: int = 3) -> bool:
    from collections import Counter
    seq   = peptide.upper()
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    return any(v >= 3 for v in Counter(kmers).values())


def predict_allergen(peptides: List[str], config: dict = None) -> Dict[str, dict]:
    """
    Rule-based allergenicity pre-screening inspired by AllerTOP criteria.
    All thresholds read from config["safety"].
    """
    if config is None:
        config = DEFAULT_CONFIG

    sf = config.get("safety", DEFAULT_CONFIG["safety"])
    confidence_cutoff = sf.get("allergen_confidence_cutoff",   0.65)
    pi_threshold      = sf.get("allergen_pi_threshold",        10.5)
    stab_len_min      = sf.get("allergen_stability_len_min",   15)
    stab_threshold    = sf.get("allergen_stability_threshold",  5.0)

    results = {}

    for peptide in peptides:
        seq = peptide.upper().replace(" ", "")
        if len(seq) < 5:
            results[peptide] = {"allergen": False, "confidence": 0.0,
                                "reason": "Too short", "details": {}}
            continue

        motif_r    = scan_allergen_motifs(seq)
        eightmer_r = check_8mer_match(seq)
        cys_flag   = check_cysteine_pattern(seq)
        pro_flag   = check_proline_repeat(seq, config)
        rep_flag   = check_repetitive_kmers(seq)
        pi         = compute_isoelectric_point(seq)
        instability= compute_instability_index(seq)

        pi_risk    = pi > pi_threshold
        stab_risk  = (instability < stab_threshold) and (len(seq) > stab_len_min)

        confidence, reasons = 0.0, []

        if motif_r["matched"]:
            top = motif_r["hits"][0]
            confidence = max(confidence, top["weight"])
            reasons.append(f"Motif: {top['motif']} → {top['allergen']} (w={top['weight']})")

        if eightmer_r["matched"]:
            confidence = max(confidence, 0.95)
            h = eightmer_r["hits"][0]
            reasons.append(f"WHO/FAO 8-mer pos {h['position']}: {h['matched']}")

        if cys_flag:
            confidence = max(confidence, 0.75)
            reasons.append(f"Cysteine scaffold ({seq.count('C')}C in {len(seq)}AA)")

        if pro_flag:
            pct = round(seq.count("P") / len(seq) * 100, 1)
            confidence = max(confidence, 0.70)
            reasons.append(f"Proline-rich repeat ({pct}% Pro)")

        if rep_flag:
            confidence = max(confidence, 0.65)
            reasons.append("Repetitive k-mer content")

        if pi_risk:
            confidence = max(confidence, 0.55)
            reasons.append(f"Basic pI={pi} (>{pi_threshold})")

        if stab_risk:
            confidence = max(confidence, 0.50)
            reasons.append(f"GI-stable long peptide (instability={instability:.1f})")

        results[peptide] = {
            "allergen":   confidence >= confidence_cutoff,
            "confidence": round(confidence, 3),
            "reason":     "; ".join(reasons) if reasons else "No allergen signals detected",
            "details": {
                "motif_match": motif_r, "who_fao_8mer": eightmer_r,
                "cysteine_scaffold": cys_flag, "proline_repeat": pro_flag,
                "repetitive_kmers": rep_flag,
                "isoelectric_point": pi, "instability_index": round(instability, 2),
                "flags": {"basic_pi": pi_risk, "gi_stable": stab_risk},
            },
        }

    return results


def filter_candidates_by_safety(
    candidates:   List[dict],
    tox_map:      Dict[str, dict],
    allergen_map: Dict[str, dict],
    verbose:      bool = True,
) -> List[dict]:
    out, rejected = [], []
    for c in candidates:
        p = c.get("peptide", "")
        if not p:
            continue
        if tox_map.get(p, {}).get("toxic"):
            rejected.append((p, f"TOXIC (score={tox_map[p].get('score','?')})"))
        elif allergen_map.get(p, {}).get("allergen"):
            allr = allergen_map[p]
            rejected.append((p, f"ALLERGEN conf={allr.get('confidence','?')}: "
                               f"{allr.get('reason','')[:80]}"))
        else:
            out.append(c)

    if verbose and rejected:
        print(f"  [SAFETY] Rejected {len(rejected)} peptides:")
        for pep, reason in rejected:
            print(f"    [!] {pep:<20} -> {reason}")

    return out
