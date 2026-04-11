# validators.py

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def validate_sequence(seq, config):
    seq = seq.upper().strip().replace("\n", "")

    if len(seq) < config["validation"]["min_length"]:
        raise ValueError("Sequence too short")

    for aa in seq:
        if aa not in VALID_AA:
            raise ValueError(f"Invalid amino acid: {aa}")

    return seq