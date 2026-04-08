# main.py

import sys
from validators import validate_sequence
from processing import get_mhci, get_mhcii, get_bcell
from bcell import extract_bcell_epitopes
from scoring import score_mhci, score_mhcii
from fusion import fuse_epitopes
from default_config import DEFAULT_CONFIG
from filter import is_redundant


def select_diverse_top(candidates, count, already_selected_peptides, overlap_threshold):
    """Selects the top N peptides that are not redundant with already selected ones."""
    selected = []
    # Candidates are already sorted by score
    for item in candidates:
        if len(selected) >= count:
            break
            
        peptide = item["peptide"]
        
        # Diversity check: not redundant with globally selected or currently selected in this batch
        if not is_redundant(peptide, already_selected_peptides + [s["peptide"] for s in selected], overlap_threshold):
            selected.append(item)
            
    return selected


def run_module1(seq, config=None):
    if config is None:
        config = DEFAULT_CONFIG

    print("🔍 Validating sequence...")
    seq = validate_sequence(seq, config)

    print("🌐 Fetching MHC-I...")
    mhci_df = get_mhci(seq, config)
    print(f"✔ MHC-I: {len(mhci_df)} hits")

    print("🌐 Fetching MHC-II...")
    mhcii_df = get_mhcii(seq, config)
    print(f"✔ MHC-II: {len(mhcii_df)} hits")

    print("🌐 Fetching B-cell...")
    bcell_df = get_bcell(seq, config)

    print("🧬 Extracting B-cell epitopes...")
    bcell_epitopes = extract_bcell_epitopes(bcell_df, config)
    print(f"✔ B-cell: {len(bcell_epitopes)} epitopes")

    print("📈 Scoring...")
    mhci_scores = score_mhci(mhci_df, config)
    mhcii_scores = score_mhcii(mhcii_df, config)

    # 1. Global Deduplication and Ranking
    # Deduplicate within categories first to get best score for each peptide
    def get_best_unique(scores):
        unique = {}
        for item in scores:
            p = item["peptide"]
            if p not in unique or item["score"] > unique[p]["score"]:
                unique[p] = item
        return sorted(unique.values(), key=lambda x: x["score"], reverse=True)

    mhci_candidates = get_best_unique(mhci_scores)
    mhcii_candidates = get_best_unique(mhcii_scores)
    
    bcell_score = config["scoring"]["bcell_score"]
    bcell_candidates = sorted([
        {"peptide": p, "type": "B-cell", "score": bcell_score} 
        for p in set(bcell_epitopes)
    ], key=lambda x: len(x["peptide"]), reverse=True)

    # 2. Diversity-Aware Selection with Quotas
    selected_all = []
    overlap_limit = config["thresholds"].get("overlap_limit", 6)
    
    # MHC-I Selection
    mhc1_selected = select_diverse_top(mhci_candidates, config["selection"]["counts"]["MHC-I"], [], overlap_limit)
    selected_all.extend(mhc1_selected)
    
    # MHC-II Selection (must not overlap with MHC-I)
    mhc2_selected = select_diverse_top(mhcii_candidates, config["selection"]["counts"]["MHC-II"], [s["peptide"] for s in selected_all], overlap_limit)
    selected_all.extend(mhc2_selected)
    
    # B-cell Selection (must not overlap with MHC-I or MHC-II)
    bcell_selected = select_diverse_top(bcell_candidates, config["selection"]["counts"]["B-cell"], [s["peptide"] for s in selected_all], overlap_limit)
    selected_all.extend(bcell_selected)

    print("🏆 Ranking and Filtering...")
    print(f"   Selected (Diverse): {len(mhc1_selected)} MHC-I, {len(mhc2_selected)} MHC-II, {len(bcell_selected)} B-cell")

    print("🔗 Building antigen...")
    antigen = fuse_epitopes(selected_all, config)

    return {
        "top_epitopes": selected_all,
        "antigen_sequence": antigen,
        "mhci_raw": mhci_df.to_dict(),
        "mhcii_raw": mhcii_df.to_dict(),
        "bcell": bcell_epitopes
    }


if __name__ == "__main__":
    if not sys.stdin.isatty():
        seq = sys.stdin.read().strip()
    else:
        seq = input("Enter protein sequence:\n").strip()
    
    if not seq:
        print("Error: No sequence provided.")
        sys.exit(1)

    try:
        result = run_module1(seq)
        print("\n=== FINAL ANTIGEN ===")
        print(result["antigen_sequence"])
    except Exception as e:
        print(f"Error occurred: {e}")
        sys.exit(1)