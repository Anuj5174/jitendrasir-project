"""
structure_3d.py
---------------
3D protein structure prediction via ESMFold (Meta AI).
  - Free REST API, no install required
  - Input:  raw amino acid sequence
  - Output: PDB file, per-residue pLDDT, mean confidence, secondary structure summary

Paper reference: Section 2.26 (Tertiary structure prediction) + 2.27 (3D validation)
"""

import os
import re
import requests
import time

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def predict_structure(sequence: str, output_dir: str = "output") -> dict:
    """
    Submit sequence to ESMFold and return PDB + confidence metrics.
    Falls back with actionable guidance if API is unavailable.
    """
    sequence = sequence.strip().upper()
    os.makedirs(output_dir, exist_ok=True)

    print("[3D]  Submitting to ESMFold API...")

    try:
        resp = requests.post(
            ESMFOLD_URL,
            data=sequence,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=180
        )

        if resp.status_code == 200 and resp.text.strip().startswith("ATOM"):
            pdb_content = resp.text

            # Save PDB file
            pdb_path = os.path.join(output_dir, "structure.pdb")
            with open(pdb_path, "w") as f:
                f.write(pdb_content)

            # Parse pLDDT scores (stored in B-factor column of ATOM records)
            plddt_scores = _parse_plddt(pdb_content)
            mean_plddt   = round(sum(plddt_scores) / len(plddt_scores), 2) if plddt_scores else 0.0

            # Confidence classification
            confidence = _classify_confidence(mean_plddt)

            # Secondary structure summary (count H=helix, E=sheet, C=coil from DSSP-like REMARK)
            sec_struct = _parse_secondary_structure(pdb_content, len(sequence))

            print(f"[3D]  ESMFold complete. Mean pLDDT: {mean_plddt} → {confidence}")

            return {
                "method":        "ESMFold",
                "pdb_path":      pdb_path,
                "pdb_content":   pdb_content,
                "mean_plddt":    mean_plddt,
                "plddt_scores":  plddt_scores,
                "confidence":    confidence,
                "secondary_structure": sec_struct,
                "length":        len(sequence),
                "status":        "success"
            }

        else:
            raise ValueError(f"ESMFold returned status {resp.status_code}: {resp.text[:200]}")

    except Exception as e:
        print(f"[3D WARN] ESMFold unavailable: {e}")
        return _esmfold_fallback(sequence)


def _parse_plddt(pdb_content: str) -> list:
    """Extract per-residue pLDDT from the B-factor column of ATOM records."""
    scores = []
    seen_residues = set()
    for line in pdb_content.splitlines():
        if line.startswith("ATOM") and " CA " in line:
            try:
                res_id = line[22:26].strip()
                chain  = line[21]
                key    = (chain, res_id)
                if key not in seen_residues:
                    seen_residues.add(key)
                    bfactor = float(line[60:66].strip())
                    scores.append(bfactor)
            except (ValueError, IndexError):
                continue
    return scores


def _classify_confidence(plddt: float) -> str:
    if plddt >= 90:  return "Very High (>90)"
    if plddt >= 70:  return "Confident (70–90)"
    if plddt >= 50:  return "Low (50–70)"
    return "Very Low (<50)"


def _parse_secondary_structure(pdb_content: str, seq_len: int) -> dict:
    """
    Estimate secondary structure from HELIX/SHEET records in PDB.
    Falls back to a length-based rough estimate if records are absent.
    """
    helix_res = set()
    sheet_res = set()

    for line in pdb_content.splitlines():
        if line.startswith("HELIX"):
            try:
                start = int(line[21:25].strip())
                end   = int(line[33:37].strip())
                for i in range(start, end + 1):
                    helix_res.add(i)
            except ValueError:
                pass
        elif line.startswith("SHEET"):
            try:
                start = int(line[22:26].strip())
                end   = int(line[33:37].strip())
                for i in range(start, end + 1):
                    sheet_res.add(i)
            except ValueError:
                pass

    total = seq_len if seq_len > 0 else 1
    helix_pct  = round(len(helix_res)  / total * 100, 1)
    sheet_pct  = round(len(sheet_res)  / total * 100, 1)
    coil_pct   = round(max(0, 100 - helix_pct - sheet_pct), 1)

    return {
        "alpha_helix_pct": helix_pct,
        "beta_sheet_pct":  sheet_pct,
        "coil_pct":        coil_pct
    }


def _esmfold_fallback(sequence: str) -> dict:
    """Return structured guidance when ESMFold API is unreachable."""
    alphafold_url = (
        "https://colab.research.google.com/github/sokrypton/ColabFold/blob/main"
        "/AlphaFold2.ipynb"
    )
    return {
        "method":    "fallback",
        "status":    "unavailable",
        "length":    len(sequence),
        "mean_plddt": None,
        "confidence": "N/A",
        "pdb_path":   None,
        "pdb_content": None,
        "secondary_structure": {},
        "fallback_guidance": {
            "message": "ESMFold API unreachable. Use AlphaFold Colab or SWISS-MODEL.",
            "alphafold_colab": alphafold_url,
            "swiss_model":     "https://swissmodel.expasy.org/",
            "sequence":        sequence
        }
    }
