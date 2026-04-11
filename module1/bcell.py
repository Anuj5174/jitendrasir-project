# bcell.py

def extract_bcell_epitopes(df, config):
    epitopes = []
    current = ""
    threshold = config["thresholds"]["bcell"]
    min_len = config["thresholds"]["min_bcell_length"]

    if df.empty:
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