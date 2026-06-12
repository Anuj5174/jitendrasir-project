"""
expression_vector.py
---------------------
In silico cloning and expression vector design (PDF §2.33 + §2.34).

Generates:
  - Cloning insert (DNA) with restriction sites: NcoI (5') and XhoI (3')
  - Forward & Reverse PCR primers with overhangs
  - Annotated sequence map
  - Vector metadata for pET-28a(+) (most common bacterial expression vector)
  - Links: C-ImmSim (§2.13), iMODS (§2.28), PROCHECK (§2.27)
"""

import re

# ── Standard molecular biology constants ─────────────────────────────────────
KOZAK        = "GCCACCATG"          # Kozak + ATG for mammalian expression
RBS          = "AAGGAG"             # Shine–Dalgarno for bacterial (E. coli)
NCOI_SITE    = "CCATGG"             # NcoI restriction site (contains ATG)
XHOI_SITE    = "CTCGAG"             # XhoI restriction site
STOP_CODON   = "TGA"                # Mammalian-optimised stop
SPACER       = "GGT"                # Glycine codon spacer between XhoI and stop
HIS_TAG_DNA  = "CATCACCATCACCATCAC" # 6×His-tag DNA (added if not in sequence)

VECTOR_DB = {
    "pET-28a(+)": {
        "host":        "E. coli (BL21-DE3)",
        "promoter":    "T7",
        "resistance":  "Kanamycin (30 µg/mL)",
        "size_kb":     5.4,
        "n_tag":       "His×6 + Thrombin site",
        "c_tag":       "His×6",
        "cloning_sites": ["NcoI", "XhoI"],
        "inducer":     "IPTG (0.5–1 mM)",
        "description": "High-copy T7-promoter vector. Best for soluble cytoplasmic expression. C-terminal His-tag auto-appended.",
        "reference":   "Novagen pET System Manual",
    },
    "pUC19": {
        "host":        "E. coli (DH5α)",
        "promoter":    "lac",
        "resistance":  "Ampicillin (100 µg/mL)",
        "size_kb":     2.7,
        "n_tag":       "None",
        "c_tag":       "None",
        "cloning_sites": ["EcoRI", "HindIII"],
        "inducer":     "IPTG (0.1 mM)",
        "description": "Cloning and propagation vector. Use pET-28a for expression.",
        "reference":   "GenBank L09137",
    },
}

CODON_TABLE = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
    'E': 'GAA', 'Q': 'CAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}


def _protein_to_dna(protein_seq: str) -> str:
    """Translate protein → DNA using E. coli preferred codons."""
    return ''.join(CODON_TABLE.get(aa.upper(), 'NNN') for aa in protein_seq)


def _reverse_complement(dna: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(dna.upper()))


def _find_restriction_sites(dna: str) -> dict:
    """Scan for common restriction sites that may interfere with cloning."""
    enzymes = {
        'NcoI':   'CCATGG',
        'XhoI':   'CTCGAG',
        'NdeI':   'CATATG',
        'BamHI':  'GGATCC',
        'EcoRI':  'GAATTC',
        'HindIII':'AAGCTT',
    }
    found = {}
    for name, site in enzymes.items():
        positions = [m.start() + 1 for m in re.finditer(site, dna)]
        if positions:
            found[name] = positions
    return found


def design_expression_vector(
    protein_seq: str,
    codon_optimised_dna: str = None,
    vector: str = "pET-28a(+)",
    primer_tm_target: int = 60
) -> dict:
    """
    Full in silico cloning pipeline.

    Args:
        protein_seq:         Vaccine construct amino acid sequence.
        codon_optimised_dna: If provided, used as the insert CDS. 
                             Otherwise reverse-translated from protein.
        vector:              Target expression vector name.
        primer_tm_target:    Desired primer melting temperature (°C).

    Returns:
        Dictionary with full cloning metadata, insert, primers, and submission links.
    """
    protein_seq = protein_seq.strip().upper().replace('*', '')

    # ── Step 1: CDS ──────────────────────────────────────────────────────────
    if codon_optimised_dna and len(codon_optimised_dna) >= len(protein_seq) * 3:
        cds = codon_optimised_dna.upper().replace('U', 'T')
        # Strip any existing stop codon
        if cds.endswith(('TGA', 'TAA', 'TAG')):
            cds = cds[:-3]
    else:
        cds = _protein_to_dna(protein_seq)

    # Ensure starts with ATG
    if not cds.startswith('ATG'):
        cds = 'ATG' + cds

    # ── Step 2: Check for internal restriction sites ──────────────────────────
    internal_sites = _find_restriction_sites(cds)

    # ── Step 3: Build complete insert ─────────────────────────────────────────
    # Structure: [NcoI]--[CDS (no stop)]--[Spacer]--[XhoI]--[STOP]
    insert = NCOI_SITE[:-3] + cds + SPACER + XHOI_SITE + STOP_CODON
    # NcoI = CCATGG → the last ATG IS the start codon, so we avoid doubling ATG:
    # Full: CC + ATG + [rest_of_cds_without_ATG] + GGT + CTCGAG + TGA
    insert_clean = "CC" + cds + SPACER + XHOI_SITE + STOP_CODON

    full_length_bp = len(insert_clean)

    # ── Step 4: PCR primers ───────────────────────────────────────────────────
    # Forward primer: NcoI overhang + first 20 nt of CDS
    fwd_overhang = "CATG" + NCOI_SITE  # 4nt clamp + NcoI
    fwd_binding  = cds[:20]
    fwd_primer   = fwd_overhang + fwd_binding

    # Reverse primer: 4nt clamp + XhoI + reverse complement of last 20 nt
    rev_overhang = "CATG" + XHOI_SITE
    rev_binding  = _reverse_complement(cds[-20:] + SPACER)
    rev_primer   = rev_overhang + rev_binding

    # Estimated Tm (Wallace rule for the binding portion only)
    def _tm(seq):
        s = seq.upper()
        return 2 * (s.count('A') + s.count('T')) + 4 * (s.count('G') + s.count('C'))

    fwd_tm = _tm(fwd_binding)
    rev_tm = _tm(rev_binding[:20])

    # ── Step 5: Vector metadata ───────────────────────────────────────────────
    vec_info = VECTOR_DB.get(vector, VECTOR_DB["pET-28a(+)"])

    # ── Step 6: Annotated sequence map ────────────────────────────────────────
    annotation_map = [
        {"label": "NcoI Site",      "start": 1,                    "end": 6,                       "sequence": "CCATGG"},
        {"label": "Start Codon",    "start": 3,                    "end": 5,                       "sequence": "ATG"},
        {"label": "CDS",            "start": 3,                    "end": 2 + len(cds),            "sequence": f"{cds[:10]}...{cds[-10:]}"},
        {"label": "Spacer",         "start": 3 + len(cds),        "end": 5 + len(cds),            "sequence": SPACER},
        {"label": "XhoI Site",      "start": 6 + len(cds),        "end": 11 + len(cds),           "sequence": "CTCGAG"},
        {"label": "Stop Codon",     "start": 12 + len(cds),       "end": 14 + len(cds),           "sequence": STOP_CODON},
    ]

    # ── Step 7: External tool links ───────────────────────────────────────────
    submission_links = {
        "C-ImmSim (Immune Simulation §2.13)": {
            "url":         "https://kraken.iac.rm.cnr.it/C-IMMSIM/index.php",
            "description": "Paste the vaccine protein sequence. Set injection schedule: Day 1, 84, 168.",
            "input":       protein_seq,
        },
        "iMODS (Molecular Dynamics §2.28)": {
            "url":         "https://imods.iqfr.csic.es/",
            "description": "Upload the downloaded PDB file to analyze RMSD, RMSF, deformability, and B-factors.",
        },
        "PROCHECK (3D Validation §2.27)": {
            "url":         "https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/",
            "description": "Upload PDB for Ramachandran plot and stereochemical quality validation.",
        },
        "ProSA (Z-score §2.27)": {
            "url":         "https://prosa.services.came.sbg.ac.at/prosa.php",
            "description": "Upload PDB to get a Z-score of structural quality.",
        },
        "HDOCK (Molecular Docking §2.18)": {
            "url":         "https://hdock.phys.hust.edu.cn/",
            "description": "Upload vaccine PDB + receptor PDB (TLR4: 4G8A, TLR2: 2Z7X, HLA: 1DUZ).",
        },
        "SnapGene Viewer (Vector Map)": {
            "url":         "https://www.snapgene.com/snapgene-viewer",
            "description": "Import pET-28a(+) from SnapGene library and insert your sequence between NcoI/XhoI.",
        },
    }

    return {
        "status":          "success",
        "vector":          vector,
        "vector_info":     vec_info,
        "protein_length":  len(protein_seq),
        "cds_length_bp":   len(cds),
        "insert_length_bp": full_length_bp,
        "cds_preview":     cds[:60] + "..." if len(cds) > 60 else cds,
        "full_insert":     insert_clean,
        "internal_restriction_sites": internal_sites,
        "primers": {
            "forward": {
                "sequence": fwd_primer,
                "length":   len(fwd_primer),
                "tm_approx": fwd_tm,
                "annotation": f"NcoI overhang + first 20nt of CDS",
            },
            "reverse": {
                "sequence": rev_primer,
                "length":   len(rev_primer),
                "tm_approx": rev_tm,
                "annotation": f"XhoI overhang + RC of last 20nt + spacer",
            },
        },
        "annotation_map":  annotation_map,
        "submission_links": submission_links,
        "cloning_protocol": [
            f"1. PCR-amplify CDS using Forward + Reverse primers above.",
            f"2. Digest PCR product with NcoI + XhoI (37°C, 60 min).",
            f"3. Digest {vector} vector with NcoI + XhoI.",
            f"4. Ligate insert into vector (T4 ligase, 16°C overnight).",
            f"5. Transform into E. coli DH5α, select on {vec_info['resistance']}.",
            f"6. Verify by colony PCR and Sanger sequencing.",
            f"7. Transfer verified plasmid to {vec_info['host']} for expression.",
            f"8. Induce with {vec_info['inducer']} at OD600 ≈ 0.6.",
        ],
    }


if __name__ == "__main__":
    import json
    test_protein = "MAKLSTDELLMAARQNLVKTPRAATVLSAPQATLVAPQATLVAPQATLVAP"
    result = design_expression_vector(test_protein)
    print(json.dumps({k: v for k, v in result.items() if k != 'full_insert'}, indent=2))
    print(f"\nInsert (first 80 bp): {result['full_insert'][:80]}...")
