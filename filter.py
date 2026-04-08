# filter.py

def filter_mhci(df, config):
    if "percentile_rank" in df.columns:
        return df[df["percentile_rank"] <= config["thresholds"]["percentile_rank"]]
    elif "rank" in df.columns:
        return df[df["rank"] <= config["thresholds"]["percentile_rank"]]
    elif "ic50" in df.columns:
        return df[df["ic50"] < config["thresholds"]["ic50"]]
    return df


def filter_mhcii(df, config):
    if "percentile_rank" in df.columns:
        return df[df["percentile_rank"] <= config["thresholds"]["percentile_rank"]]
    elif "rank" in df.columns:
        return df[df["rank"] <= config["thresholds"]["percentile_rank"]]
    elif "ic50" in df.columns:
        return df[df["ic50"] < config["thresholds"]["ic50"]]
    return df


def get_longest_common_substring(s1, s2):
    m, n = len(s1), len(s2)
    lcs = 0
    for i in range(m):
        for j in range(n):
            k = 0
            while (i + k < m) and (j + k < n) and s1[i+k] == s2[j+k]:
                k += 1
            lcs = max(lcs, k)
    return lcs

def is_redundant(new_peptide, existing_peptides, overlap_threshold):
    """
    Returns True if the new_peptide is redundant with any already selected peptide.
    Redundancy defined as:
    1. Exact match
    2. new_peptide is a substring of an existing peptide (e.g., 'ABC' in 'ABCD')
    3. An existing peptide is a substring of new_peptide (e.g., 'ABCD' covers 'ABC')
    4. Biological overlap (sharing a motif of length >= overlap_threshold)
    """
    for existing in existing_peptides:
        if new_peptide in existing or existing in new_peptide:
            return True
        if get_longest_common_substring(new_peptide, existing) >= overlap_threshold:
            return True
    return False