"""
Standalone runner — test adjuvant auto-selection with current CSV epitopes.
Uses adjuvant_selector.py (same logic as main pipeline).
"""

import sys, os, json, csv
sys.path.insert(0, os.path.dirname(__file__))

from adjuvant_selector import auto_select_adjuvant, print_adjuvant_comparison, CANDIDATE_ADJUVANTS
from default_config import DEFAULT_CONFIG


def load_epitopes_from_csv(output_dir):
    """Dynamically read epitopes from the pipeline-generated CSV files."""
    epitopes = []

    # CTL
    ctl_path = os.path.join(output_dir, "epitope_report_ctl.csv")
    if os.path.exists(ctl_path):
        with open(ctl_path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                pep = row.get("Peptide", "").strip()
                if pep:
                    epitopes.append({"peptide": pep, "type": "MHC-I"})
        print(f"  [+] Loaded {sum(1 for e in epitopes if e['type']=='MHC-I')} CTL epitopes")

    # HTL
    htl_path = os.path.join(output_dir, "epitope_report_htl.csv")
    if os.path.exists(htl_path):
        with open(htl_path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                pep = row.get("Peptide", "").strip()
                if pep:
                    epitopes.append({"peptide": pep, "type": "MHC-II"})
        print(f"  [+] Loaded {sum(1 for e in epitopes if e['type']=='MHC-II')} HTL epitopes")

    # B-cell
    bcell_path = os.path.join(output_dir, "epitope_report_bcell.csv")
    if os.path.exists(bcell_path):
        with open(bcell_path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                pep = row.get("Peptide", "").strip()
                if pep:
                    epitopes.append({"peptide": pep, "type": "B-cell"})
        print(f"  [+] Loaded {sum(1 for e in epitopes if e['type']=='B-cell')} B-cell epitopes")

    if not epitopes:
        print("\n  [ERROR] No epitopes found! Run the pipeline first.")
        sys.exit(1)

    return epitopes


if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(__file__), "output")
    print("\n  Loading epitopes from pipeline CSV output...\n")
    epitopes = load_epitopes_from_csv(output_dir)

    print(f"\n  Total: {len(epitopes)} epitopes")
    print("  Running adjuvant auto-selection...\n")

    result = auto_select_adjuvant(epitopes, DEFAULT_CONFIG)

    print("=" * 85)
    print("  ADJUVANT AUTO-SELECTION REPORT")
    print("=" * 85)
    print_adjuvant_comparison(result)

    # Print winner construct
    w = result["all_results"][result["winner_name"]]
    print("=" * 85)
    print(f"  FINAL VACCINE CONSTRUCT ({result['winner_name']})")
    print("=" * 85)
    print(f"\n  Length: {w['sequence_length']} aa  |  MW: {w['molecular_weight_kda']} kDa")
    print(f"  Antigenicity: {w['antigenicity_score']:.4f}  |  II: {w['instability_index']:.2f} ({w['stability']})")
    print(f"  GRAVY: {w['gravy']:.3f} ({w['solubility']})\n")
    seq = w["construct"]
    for j in range(0, len(seq), 60):
        print(f"  {seq[j:j+60]}")

    # Save
    out_path = os.path.join(output_dir, "vaccine_construct_comparison.json")
    export = {}
    for k, v in result["all_results"].items():
        export[k] = {kk: vv for kk, vv in v.items() if kk != "construct"}
    with open(out_path, "w") as f:
        json.dump({
            "ranking": [{"rank": i+1, "adjuvant": n, "score": s} for i, (n, s) in enumerate(result["ranking"])],
            "details": export,
            "weights": result["weights"],
            "winner": result["winner_name"],
        }, f, indent=2)
    print(f"\n  Saved -> {out_path}")

    fasta_path = os.path.join(output_dir, "final_vaccine_construct.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">VaccineConstruct|{result['winner_name']}|{w['sequence_length']}aa|{w['molecular_weight_kda']}kDa\n")
        for j in range(0, len(seq), 70):
            f.write(seq[j:j+70] + "\n")
    print(f"  Saved -> {fasta_path}")
