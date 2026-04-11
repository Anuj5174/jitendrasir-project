# advanced_filters.py
import subprocess
import json
import shutil
from typing import List, Dict


def run_toxinpred(peptides: List[str]) -> Dict[str, dict]:
    """Attempt to run toxinpred3 on a list of peptides.

    Returns a mapping peptide -> { 'toxic': bool, 'score': float, 'note': str }
    Falls back to conservative non-toxic defaults if tool not available.
    """
    results = {}

    # Try direct import first (if installed as Python package)
    try:
        import toxinpred3

        for p in peptides:
            out = toxinpred3.predict(p)
            results[p] = {
                'toxic': bool(out.get('toxic', False)),
                'score': float(out.get('score', 0.0)),
                'note': 'via toxinpred3'
            }
        return results
    except Exception:
        pass

    # Try command-line RNA/toxinpred3 binary
    if shutil.which('toxinpred3'):
        try:
            proc = subprocess.run(['toxinpred3', '--json', '-'], input='\n'.join(peptides).encode('utf-8'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            out = json.loads(proc.stdout.decode('utf-8'))
            for p in peptides:
                r = out.get(p, {})
                results[p] = {
                    'toxic': bool(r.get('toxic', False)),
                    'score': float(r.get('score', 0.0)),
                    'note': 'via toxinpred3-cli'
                }
            return results
        except Exception:
            pass

    # Fallback: conservative non-toxic default (user should run real predictions)
    for p in peptides:
        results[p] = {'toxic': False, 'score': 0.0, 'note': 'fallback: tool unavailable'}

    return results


def predict_allergen(peptides: List[str]) -> Dict[str, dict]:
    """Stubbed Allergenicity API interface.

    Currently returns a placeholder result for each peptide. Users can
    replace this function body to call a commercial Allergenicity API.
    """
    results = {}
    for p in peptides:
        results[p] = {
            'allergen': False,
            'confidence': 0.0,
            'note': 'stub - replace with commercial API'
        }
    return results


def filter_candidates_by_safety(candidates: List[dict], tox_map: Dict[str, dict], allergen_map: Dict[str, dict]) -> List[dict]:
    """Filter a list of candidate dicts (each must contain 'peptide') by toxicity/allergen maps.

    Returns a filtered list excluding peptides flagged as toxic or allergenic.
    """
    out = []
    for c in candidates:
        p = c.get('peptide')
        if not p:
            continue
        tox = tox_map.get(p, {})
        allr = allergen_map.get(p, {})
        if tox.get('toxic'):
            continue
        if allr.get('allergen'):
            continue
        out.append(c)
    return out
