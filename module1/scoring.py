# scoring.py

import math


def compute_score_from_percentile(percentile_rank, max_score_fallback=10.0):
    """Lower percentile_rank = stronger binder = higher score."""
    if percentile_rank <= 0:
        return max_score_fallback
    return -math.log(percentile_rank / 100.0)


def compute_score_from_ic50(ic50, max_score_fallback=10.0):
    """Lower ic50 = stronger binder = higher score.
    Uses -log(IC50) as requested."""
    if ic50 <= 0:
        return max_score_fallback
    return -math.log(ic50)


def _get_score_col(df):
    """Determine which scoring column is available."""
    if "ic50" in df.columns:
        return "ic50"
    elif "percentile_rank" in df.columns:
        return "percentile_rank"
    elif "rank" in df.columns:
        return "rank"
    elif "score" in df.columns:
        return "score"
    return None


def _compute_peptide_score(group, score_col, config):
    """Compute aggregate score for a peptide group.
    Improved: score = -log(IC50) + allele_count (when IC50 available).
    Falls back to percentile-based scoring otherwise."""
    max_score_fallback = config["scoring"].get("max_score_fallback", 10.0)
    if score_col == "ic50":
        # Primary formula: score = sum(-log(IC50)) + allele_count
        binding_score = sum(compute_score_from_ic50(v, max_score_fallback) for v in group[score_col])
        allele_count = len(group["allele"].unique())
        return binding_score + allele_count
    elif score_col in ("percentile_rank", "rank"):
        binding_score = sum(compute_score_from_percentile(v, max_score_fallback) for v in group[score_col])
        allele_count = len(group["allele"].unique())
        return binding_score + allele_count
    elif score_col == "score":
        binding_score = sum(compute_score_from_percentile(v * 100, max_score_fallback) for v in group[score_col])
        allele_count = len(group["allele"].unique())
        return binding_score + allele_count
    return 0


def score_mhci(df, config):
    if df.empty:
        return []
    score_col = _get_score_col(df)
    if not score_col:
        return []
    grouped = df.groupby("peptide")
    results = []
    for peptide, group in grouped:
        score = _compute_peptide_score(group, score_col, config)
        results.append({
            "peptide": peptide,
            "type": "MHC-I",
            "score": score
        })
    return results


def score_mhcii(df, config):
    if df.empty:
        return []
    score_col = _get_score_col(df)
    if not score_col:
        return []
    grouped = df.groupby("peptide")
    results = []
    for peptide, group in grouped:
        score = _compute_peptide_score(group, score_col, config)
        results.append({
            "peptide": peptide,
            "type": "MHC-II",
            "score": score
        })
    return results