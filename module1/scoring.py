# scoring.py
import math


def compute_score_from_percentile(percentile_rank, max_score_fallback=10.0):
    """Lower percentile_rank = stronger binder = higher score."""
    if percentile_rank <= 0:
        return max_score_fallback
    return -math.log(percentile_rank / 100.0)


def compute_score_from_ic50(ic50, config):
    """
    Score a peptide from its IC50 value.
    Formula: log(ic50_normalizer / IC50) + tier_bonus + allele_breadth
    All constants come from config["scoring"].
    """
    sc = config["scoring"]
    max_fallback  = sc.get("max_score_fallback",        10.0)
    normalizer    = sc.get("ic50_normalizer",         50000.0)
    strong_thr    = sc.get("ic50_strong_threshold",      50.0)
    moderate_thr  = sc.get("ic50_moderate_threshold",   500.0)
    strong_bonus  = sc.get("ic50_strong_bonus",           5.0)
    moderate_bonus= sc.get("ic50_moderate_bonus",         2.0)

    if ic50 <= 0:
        return max_fallback

    base = math.log(normalizer / max(1.0, ic50))

    tier_bonus = 0.0
    if ic50 < strong_thr:
        tier_bonus = strong_bonus
    elif ic50 < moderate_thr:
        tier_bonus = moderate_bonus

    return base + tier_bonus


def _get_score_col(df):
    """Determine which scoring column is available in the DataFrame."""
    for col in ("ic50", "percentile_rank", "rank", "score"):
        if col in df.columns:
            return col
    return None


def _compute_peptide_score(group, score_col, config):
    """
    Compute a composite score for a peptide group (one peptide, multiple alleles).

    Uses the BEST (minimum) IC50 across alleles as the binding signal,
    then adds an allele-breadth bonus so peptides binding many alleles
    rank higher than narrow binders with the same best IC50.
    """
    sc = config["scoring"]
    max_fallback    = sc.get("max_score_fallback",    10.0)
    breadth_weight  = sc.get("allele_breadth_weight",  1.0)

    allele_count  = len(group["allele"].unique())
    breadth_bonus = breadth_weight * math.log1p(allele_count)

    if score_col == "ic50":
        best_ic50     = group[score_col].min()
        binding_score = compute_score_from_ic50(best_ic50, config)
        return binding_score + breadth_bonus

    elif score_col in ("percentile_rank", "rank"):
        best_rank     = group[score_col].min()
        binding_score = compute_score_from_percentile(best_rank, max_fallback)
        return binding_score + breadth_bonus

    elif score_col == "score":
        best_val      = group[score_col].min()
        binding_score = compute_score_from_percentile(best_val * 100, max_fallback)
        return binding_score + breadth_bonus

    return 0.0


def score_mhci(df, config):
    if df.empty:
        return []
    score_col = _get_score_col(df)
    if not score_col:
        return []
    results = []
    for peptide, group in df.groupby("peptide"):
        score = _compute_peptide_score(group, score_col, config)
        results.append({"peptide": peptide, "type": "MHC-I", "score": score})
    return results


def score_mhcii(df, config):
    if df.empty:
        return []
    score_col = _get_score_col(df)
    if not score_col:
        return []
    results = []
    for peptide, group in df.groupby("peptide"):
        score = _compute_peptide_score(group, score_col, config)
        results.append({"peptide": peptide, "type": "MHC-II", "score": score})
    return results