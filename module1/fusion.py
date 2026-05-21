# fusion.py

def merge_junction(seq1, seq2):
    """
    Merge seq1 and seq2 smoothly. If seq1 ends with a sequence that seq2 starts with,
    overlap them so there is no stacking (e.g., '...TAG' + 'KK' -> '...TAGKK', 
    but '...TAGK' + 'KK' -> '...TAGKK').
    """
    max_overlap = min(len(seq1), len(seq2))
    for i in range(max_overlap, 0, -1):
        if seq1.endswith(seq2[:i]):
            return seq1 + seq2[i:]
    return seq1 + seq2

def fuse_epitopes(epitopes, config):
    """Build the complete vaccine construct following PDF-described architecture:
    
    [Adjuvant]-EAAAK-[PADRE]-EAAAK-[MHC-I epitopes joined by AAY]-KK-
    [MHC-II epitopes joined by GPGPG]-KK-[B-cell epitopes joined by KK]-HHHHHH
    
    This matches the construct design from:
    - Martinelli (2022): Adjuvant + EAAAK + PADRE + linker-separated epitopes
    - Mortazavi et al. (2024): Adjuvant at N-terminal via EAAAK
    - Shabbir et al. (2025) Nipah paper: 50S ribosomal + EAAAK + epitopes + 6xHis
    """
    if not epitopes:
        return ""

    construct_cfg = config.get("construct", {})
    linkers = construct_cfg.get("linkers", {})
    
    separator = linkers.get("B-cell", "KK")
    rigid_linker = linkers.get("rigid", "EAAAK")
    adjuvant_seq = construct_cfg.get("adjuvant_sequence", "")
    padre_seq = construct_cfg.get("padre_sequence", "")
    his_tag = construct_cfg.get("his_tag", "")

    # --- 1. Start with Adjuvant ---
    sequence = ""
    if adjuvant_seq:
        sequence = adjuvant_seq
        sequence = merge_junction(sequence, rigid_linker)

    # --- 2. Add PADRE universal T-helper peptide ---
    if padre_seq:
        sequence = merge_junction(sequence, padre_seq)
        sequence = merge_junction(sequence, rigid_linker)

    # --- 3. Group epitopes by type for ordered assembly ---
    groups = {"MHC-I": [], "MHC-II": [], "B-cell": []}
    for ep in epitopes:
        ep_type = ep.get("type", "MHC-I")
        if ep_type in groups:
            groups[ep_type].append(ep)
        else:
            groups["MHC-I"].append(ep)

    # --- 4. Fuse each group with its type-specific linker ---
    group_sequences = []
    for group_type in ["MHC-I", "MHC-II", "B-cell"]:
        group = groups[group_type]
        if not group:
            continue
        group_linker = linkers.get(group_type, "")
        group_seq = group[0]["peptide"]
        for i in range(1, len(group)):
            group_seq = merge_junction(group_seq, group_linker)
            group_seq = merge_junction(group_seq, group[i]["peptide"])
        group_sequences.append(group_seq)

    # --- 5. Join groups with separator (KK) ---
    for i, gs in enumerate(group_sequences):
        if i == 0:
            sequence = merge_junction(sequence, gs)
        else:
            sequence = merge_junction(sequence, separator)
            sequence = merge_junction(sequence, gs)

    # --- 6. Add 6x Histidine tag at C-terminus ---
    if his_tag:
        sequence = merge_junction(sequence, his_tag)

    return sequence