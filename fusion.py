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
    if not epitopes:
        return ""

    linkers = config["fusion"]["linkers"]
    sequence = epitopes[0]["peptide"]

    for i in range(1, len(epitopes)):
        current_type = epitopes[i]["type"]
        linker = linkers.get(current_type, "")
        next_peptide = epitopes[i]["peptide"]
        
        # Merge the sequence with the linker, then with the next peptide
        sequence = merge_junction(sequence, linker)
        sequence = merge_junction(sequence, next_peptide)

    return sequence