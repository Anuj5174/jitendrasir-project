# bcell.py
import pandas as pd


def extract_bcell_epitopes(df, config):
    epitopes = []
    current = ""
    threshold = config["thresholds"]["bcell"]
    min_len = config["thresholds"]["min_bcell_length"]

    if df.empty:
        print("  [INFO] No valid B-cell epitopes found under current thresholds (Linear BepiPred).")
        print("  [INFO] Recommended: Use ElliPro (structural B-cell prediction) for conformational epitopes ")
        print("         using solvent accessibility and antibody exposure scoring, as per standard protocols.")
        return []

    score_col = "Score" if "Score" in df.columns else "score" if "score" in df.columns else None
    res_col = "Residue" if "Residue" in df.columns else "residue" if "residue" in df.columns else None

    if not score_col or not res_col:
        return []

    for _, row in df.iterrows():
        if row[score_col] >= threshold:
            current += row[res_col]
        else:
            if len(current) >= min_len:
                epitopes.append(current)
            current = ""

    if len(current) >= min_len:
        epitopes.append(current)

    return epitopes


def extract_bcell_epitopes_relaxed(df, config, min_count=1):
    """
    Relaxed version of epitope extraction for cases where no epitopes meet
    the primary biological threshold. It lowers the bar to find 'candidate' regions.
    """
    if df.empty:
        print("  [INFO] No valid B-cell epitopes found even under relaxed thresholds.")
        print("  [INFO] Recommended: Use ensemble B-cell prediction or structural tools like ElliPro.")
        return []

    score_col = "Score" if "Score" in df.columns else "score" if "score" in df.columns else None
    res_col = "Residue" if "Residue" in df.columns else "residue" if "residue" in df.columns else None
    if not score_col or not res_col:
        return []

    # Strategy: Find top scoring continuous region of length >= min_len
    min_len = config["thresholds"]["min_bcell_length"]
    
    # Try with a much lower threshold first (half of original)
    relaxed_threshold = config["thresholds"]["bcell"] * 0.5
    
    epitopes = []
    current = ""
    for _, row in df.iterrows():
        if row[score_col] >= relaxed_threshold:
            current += row[res_col]
        else:
            if len(current) >= min_len:
                epitopes.append(current)
            current = ""
    if len(current) >= min_len:
        epitopes.append(current)

    # No synthetic fabrication of epitopes (Strict adherence to research protocols)
    if not epitopes:
        print("  [WARN] Adaptive thresholding failed to yield B-cell epitopes. Returning empty.")

    return epitopes[:min_count]