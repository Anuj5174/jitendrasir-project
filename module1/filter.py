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


def get_max_contiguous_match(s1, s2):
    """Return length of the longest contiguous exact match between s1 and s2."""
    m, n = len(s1), len(s2)
    best = 0
    for i in range(m):
        for j in range(n):
            k = 0
            while (i + k < m) and (j + k < n) and s1[i + k] == s2[j + k]:
                k += 1
            if k > best:
                best = k
    return best


def kmer_overlap_fraction(s1, s2, k=3):
    """Compute fraction of k-mers in s1 that also appear in s2.

    This gives a lightweight, sliding-window overlap estimate useful
    for biological motif similarity checks.
    """
    if len(s1) < k or len(s2) < k:
        return 0.0
    s2_kmers = set(s2[i:i + k] for i in range(len(s2) - k + 1))
    s1_kmers = [s1[i:i + k] for i in range(len(s1) - k + 1)]
    matches = sum(1 for x in s1_kmers if x in s2_kmers)
    return matches / max(1, len(s1_kmers))


def get_overlap_score(s1, s2, window=6):
    """
    Returns the maximum overlap length between two sequences using a sliding window.
    This identifies if s1 and s2 share a common motif of at least 'window' size.
    """
    m, n = len(s1), len(s2)
    if m < window or n < window:
        return 0
    
    # Simple k-mer set intersection for fast motif check
    def get_kmers(s, k):
        return {s[i:i+k] for i in range(len(s) - k + 1)}
    
    kmers_1 = get_kmers(s1, window)
    kmers_2 = get_kmers(s2, window)
    
    common = kmers_1.intersection(kmers_2)
    if not common:
        return 0
        
    # If we have common k-mers, find the longest exact match
    # (Existing LCS logic but we already know there's at least 'window' match)
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
    Uses a fast sliding-window overlap check.
    """
    for existing in existing_peptides:
        if new_peptide == existing or new_peptide in existing or existing in new_peptide:
            return True

        # contiguous match check
        max_contig = get_max_contiguous_match(new_peptide, existing)
        if max_contig >= overlap_threshold:
            return True

        # k-mer sliding-window overlap fraction (sensitive to motif sharing)
        frac = kmer_overlap_fraction(new_peptide, existing, k=3)
        if frac >= 0.6:
            return True

    return False