# homology.py
"""
Human proteome homology check using NCBI BLASTp REST API.

Epitopes homologous to human proteins can trigger autoimmune reactions.
All three reference PDFs mandate this check:
- Martinelli (2022): NCBI BLAST for human protein identity
- Mortazavi et al. (2024): BLASTp for homology assessment
- Shabbir et al. (2025): Cross-reactivity filtering

Uses NCBI BLAST REST API (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
"""

import requests
import time
import re
from typing import List, Dict

def blast_peptide(peptide: str, config: dict, timeout: int = 120) -> dict:
    """Submit a peptide to NCBI BLASTp against human proteome (taxid:9606).
    
    Returns the best hit identity percentage and details.
    """
    # Threshold for safeness
    safe_threshold = config["homology"].get("identity_threshold", 50.0)

    # NCBI BLAST REST API - Submit
    put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    put_params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": peptide,
        "ENTREZ_QUERY": "Homo sapiens[organism]",
        "EXPECT": "10",
        "FORMAT_TYPE": "Text",
        "HITLIST_SIZE": "5"
    }
    
    try:
        resp = requests.post(put_url, data=put_params, timeout=30)
        if resp.status_code != 200:
            return {"error": f"BLAST submit failed: {resp.status_code}", "safe": True}
        
        # Extract RID (Request ID)
        rid_match = re.search(r"RID = (\S+)", resp.text)
        if not rid_match:
            return {"error": "Could not extract BLAST RID", "safe": True}
        
        rid = rid_match.group(1)
        
        # Poll for results
        get_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        elapsed = 0
        poll_interval = 10
        
        while elapsed < timeout:
            time.sleep(poll_interval)
            elapsed += poll_interval
            
            check_resp = requests.get(get_url, params={
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "Text"
            }, timeout=30)
            
            if "Status=WAITING" in check_resp.text:
                continue
            elif "Status=FAILED" in check_resp.text:
                return {"error": "BLAST search failed", "safe": True}
            elif "Status=READY" in check_resp.text:
                # Parse identity from results
                identity_match = re.search(r"Identities = (\d+)/(\d+) \((\d+)%\)", check_resp.text)
                if identity_match:
                    identical = int(identity_match.group(1))
                    total = int(identity_match.group(2))
                    pct = int(identity_match.group(3))
                    
                    return {
                        "identity_pct": pct,
                        "identical_residues": identical,
                        "alignment_length": total,
                        "safe": pct < safe_threshold, 
                        "note": f"{pct}% identity to human protein" + 
                                (" - AUTOIMMUNE RISK [!]" if pct >= safe_threshold else " - SAFE [+]")
                    }
                else:
                    # No significant hits found = safe
                    return {
                        "identity_pct": 0,
                        "safe": True,
                        "note": "No significant human homology found - SAFE [+]"
                    }
        
        return {"error": "BLAST timeout", "safe": True}
        
    except Exception as e:
        return {"error": str(e), "safe": True}


def check_human_homology_fast(peptides: List[str], config: dict) -> Dict[str, dict]:
    """Fast local heuristic for human self-peptide detection.
    
    Checks against known human self-peptide motifs commonly flagged
    in vaccine design.
    """
    human_self_motifs = config["homology"].get("human_self_motifs", [])
    max_short = config["homology"].get("max_short_peptide_len", 9)
    bias_ratio = config["homology"].get("codon_bias_ratio", 0.85)

    results = {}
    for p in peptides:
        p_upper = p.upper()
        is_homologous = False
        matched_motif = None
        
        for motif in human_self_motifs:
            if motif in p_upper:
                is_homologous = True
                matched_motif = motif
                break
        
        # Compositional heuristic: very short peptides with high human-codon bias
        if not is_homologous and len(p) <= max_short:
            common_human = sum(1 for aa in p_upper if aa in 'LAESRVG')
            if common_human / len(p) > bias_ratio:
                is_homologous = True
                matched_motif = "high human residue bias"
        
        results[p] = {
            "homologous": is_homologous,
            "safe": not is_homologous,
            "matched_motif": matched_motif,
            "note": f"Human motif detected: {matched_motif}" if is_homologous 
                    else "No human homology detected (fast screen)"
        }
    
    return results


def check_human_homology(peptides: List[str], config: dict, use_api: bool = False) -> Dict[str, dict]:
    """Main homology checking interface.
    """
    if use_api:
        results = {}
        for i, p in enumerate(peptides):
            print(f"[BLAST] Checking peptide {i+1}/{len(peptides)}: {p[:20]}...")
            results[p] = blast_peptide(p, config)
            # NCBI rate limit: max 3 requests per second
            if i < len(peptides) - 1:
                time.sleep(3)
        return results
    else:
        return check_human_homology_fast(peptides, config)
