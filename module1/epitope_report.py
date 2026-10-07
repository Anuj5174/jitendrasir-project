# epitope_report.py
"""
CSV and JSON export for selected epitopes with all multi-criteria properties.
Generates publication-ready tables for the vaccine design pipeline.
"""

import csv
import json
import os
from typing import List
from datetime import datetime


# Column definitions for the CTL (CD8+ / MHC-I) epitope report
CTL_COLUMNS = [
    "Rank", "Peptide", "Length", "Type",
    "Binding_Score", "Best_IC50_nM",
    "Antigenicity_Score", "Is_Antigenic",
    "Allergenicity", "Toxicity",
    "Immunogenicity_Score", "Immunogenicity_Grade",
    "TC50_nM", "TC50_Grade",
    "HLA_Allele_Count", "HLA_Alleles",
    "Population_Coverage_%", "Conservancy_%",
    "Composite_Rank",
    "Binding_Norm", "Antigenicity_Norm", "Allergenicity_Norm",
    "Toxicity_Norm", "Immunogenicity_Norm", "TC50_Norm", "MultiHLA_Norm",
]

# Column definitions for HTL (CD4+ / MHC-II) — same 7 criteria + 3 cytokines
HTL_COLUMNS = [
    "Rank", "Peptide", "Length", "Type",
    "Binding_Score", "Best_IC50_nM",
    "Antigenicity_Score", "Is_Antigenic",
    "Allergenicity", "Toxicity",
    "Immunogenicity_Score", "Immunogenicity_Grade",
    "TC50_nM", "TC50_Grade",
    "HLA_Allele_Count", "HLA_Alleles",
    "IFN_Gamma_Score", "IFN_Gamma_Inducer",
    "IL4_Score", "IL4_Inducer",
    "IL10_Score", "IL10_Inducer",
    "Th_Bias",
    "Population_Coverage_%", "Conservancy_%",
    "Composite_Rank",
    "Binding_Norm", "Antigenicity_Norm", "Allergenicity_Norm",
    "Toxicity_Norm", "Immunogenicity_Norm", "TC50_Norm", "MultiHLA_Norm",
    "IFN_Gamma_Norm", "IL4_Norm", "IL10_Norm",
]

# Column definitions for B-cell epitopes
BCELL_COLUMNS = [
    "Rank", "Peptide", "Length", "Type",
    "Binding_Score",
    "Antigenicity_Score", "Is_Antigenic",
    "Allergenicity", "Toxicity",
    "Population_Coverage_%", "Conservancy_%",
    "Composite_Rank",
    "Antigenicity_Norm", "Allergenicity_Norm", "Toxicity_Norm",
]

# Combined report uses the HTL superset (widest column set)
ALL_COLUMNS = list(HTL_COLUMNS)


def _ctl_row(rank: int, ep: dict) -> dict:
    """Build a single CTL row dict from an epitope dict."""
    sub = ep.get("_sub", {})
    alleles = ep.get("alleles_bound", [])
    return {
        "Rank": rank,
        "Peptide": ep.get("peptide", ""),
        "Length": len(ep.get("peptide", "")),
        "Type": "CTL (CD8+ / MHC-I)",
        "Binding_Score": round(ep.get("score", 0), 4),
        "Best_IC50_nM": ep.get("best_ic50", "N/A"),
        "Antigenicity_Score": round(ep.get("antigenicity", 0), 4),
        "Is_Antigenic": "YES" if ep.get("is_antigen") else "NO",
        "Allergenicity": "NON-ALLERGEN" if not ep.get("is_allergen") else "ALLERGEN",
        "Toxicity": "NON-TOXIC" if not ep.get("is_toxic") else "TOXIC",
        "Immunogenicity_Score": round(ep.get("immunogenicity_score", 0), 4),
        "Immunogenicity_Grade": ep.get("immunogenicity_grade", "N/A"),
        "TC50_nM": round(ep.get("tc50_nM", 0), 2),
        "TC50_Grade": ep.get("tc50_grade", "N/A"),
        "HLA_Allele_Count": ep.get("allele_count", 0),
        "HLA_Alleles": "; ".join(alleles),
        "IFN_Gamma_Score": "", "IFN_Gamma_Inducer": "",
        "IL4_Score": "", "IL4_Inducer": "",
        "IL10_Score": "", "IL10_Inducer": "",
        "Th_Bias": "",
        "Population_Coverage_%": round(ep.get("population_coverage", 0) * 100, 2),
        "Conservancy_%": round(ep.get("conservancy", 0) * 100, 2),
        "Composite_Rank": round(ep.get("composite_rank", 0), 4),
        "Binding_Norm": sub.get("binding", 0),
        "Antigenicity_Norm": sub.get("antigenicity", 0),
        "Allergenicity_Norm": sub.get("allergenicity", 0),
        "Toxicity_Norm": sub.get("toxicity", 0),
        "Immunogenicity_Norm": sub.get("immunogenicity", 0),
        "TC50_Norm": sub.get("tc50", 0),
        "MultiHLA_Norm": sub.get("multi_hla", 0),
        "IFN_Gamma_Norm": "", "IL4_Norm": "", "IL10_Norm": "",
    }


def _htl_row(rank: int, ep: dict) -> dict:
    """Build a row dict for HTL epitopes with all 10 criteria."""
    sub = ep.get("_sub", {})
    alleles = ep.get("alleles_bound", [])
    return {
        "Rank": rank,
        "Peptide": ep.get("peptide", ""),
        "Length": len(ep.get("peptide", "")),
        "Type": "HTL (CD4+ / MHC-II)",
        "Binding_Score": round(ep.get("score", 0), 4),
        "Best_IC50_nM": ep.get("best_ic50", "N/A"),
        "Antigenicity_Score": round(ep.get("antigenicity", 0), 4),
        "Is_Antigenic": "YES" if ep.get("is_antigen") else "NO",
        "Allergenicity": "NON-ALLERGEN" if not ep.get("is_allergen") else "ALLERGEN",
        "Toxicity": "NON-TOXIC" if not ep.get("is_toxic") else "TOXIC",
        "Immunogenicity_Score": round(ep.get("immunogenicity_score", 0), 4),
        "Immunogenicity_Grade": ep.get("immunogenicity_grade", "N/A"),
        "TC50_nM": round(ep.get("tc50_nM", 0), 2),
        "TC50_Grade": ep.get("tc50_grade", "N/A"),
        "HLA_Allele_Count": ep.get("allele_count", 0),
        "HLA_Alleles": "; ".join(alleles),
        "IFN_Gamma_Score": round(ep.get("ifng_score", 0), 4),
        "IFN_Gamma_Inducer": "YES" if ep.get("ifng_positive") else "NO",
        "IL4_Score": round(ep.get("il4_score", 0), 4),
        "IL4_Inducer": "YES" if ep.get("il4_positive") else "NO",
        "IL10_Score": round(ep.get("il10_score", 0), 4),
        "IL10_Inducer": "YES" if ep.get("il10_positive") else "NO",
        "Th_Bias": ep.get("th_bias", "Unknown"),
        "Population_Coverage_%": round(ep.get("population_coverage", 0) * 100, 2),
        "Conservancy_%": round(ep.get("conservancy", 0) * 100, 2),
        "Composite_Rank": round(ep.get("composite_rank", 0), 4),
        "Binding_Norm": sub.get("binding", 0),
        "Antigenicity_Norm": sub.get("antigenicity", 0),
        "Allergenicity_Norm": sub.get("allergenicity", 0),
        "Toxicity_Norm": sub.get("toxicity", 0),
        "Immunogenicity_Norm": sub.get("immunogenicity", 0),
        "TC50_Norm": sub.get("tc50", 0),
        "MultiHLA_Norm": sub.get("multi_hla", 0),
        "IFN_Gamma_Norm": sub.get("ifn_gamma", 0),
        "IL4_Norm": sub.get("il4", 0),
        "IL10_Norm": sub.get("il10", 0),
    }


def _bcell_row(rank: int, ep: dict) -> dict:
    """Build a row dict for B-cell epitopes."""
    sub = ep.get("_sub", {})
    return {
        "Rank": rank,
        "Peptide": ep.get("peptide", ""),
        "Length": len(ep.get("peptide", "")),
        "Type": "B-cell (Humoral)",
        "Binding_Score": round(ep.get("score", 0), 4),
        "Antigenicity_Score": round(ep.get("antigenicity", 0), 4),
        "Is_Antigenic": "YES" if ep.get("is_antigen") else "NO",
        "Allergenicity": "NON-ALLERGEN" if not ep.get("is_allergen") else "ALLERGEN",
        "Toxicity": "NON-TOXIC" if not ep.get("is_toxic") else "TOXIC",
        "Population_Coverage_%": round(ep.get("population_coverage", 0) * 100, 2),
        "Conservancy_%": round(ep.get("conservancy", 0) * 100, 2),
        "Composite_Rank": round(ep.get("composite_rank", 0), 4),
        "Antigenicity_Norm": sub.get("antigenicity", 0),
        "Allergenicity_Norm": sub.get("allergenicity", 0),
        "Toxicity_Norm": sub.get("toxicity", 0),
    }


def export_epitope_csv(
    ctl_epitopes: List[dict],
    htl_epitopes: List[dict],
    bcell_epitopes: List[dict],
    combined_ctl_coverage: float,
    output_dir: str = "output",
) -> dict:
    """
    Export all selected epitopes to CSV files.

    Generates:
      1. epitope_report_ctl.csv   — CD8+/MHC-I with 7 criteria
      2. epitope_report_htl.csv   — CD4+/MHC-II with 10 criteria (7 + cytokines)
      3. epitope_report_bcell.csv — B-cell epitopes
      4. epitope_report_all.csv   — Combined table (all types)

    Returns dict with file paths and summary data for the API.
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    paths = {}

    # ── 1. CTL CSV ────────────────────────────────────────────────────────────
    ctl_path = os.path.join(output_dir, "epitope_report_ctl.csv")
    ctl_rows = [_ctl_row(i, ep) for i, ep in enumerate(ctl_epitopes, 1)]
    _write_csv(ctl_path, CTL_COLUMNS, ctl_rows)
    paths["ctl_csv"] = ctl_path

    # ── 2. HTL CSV ────────────────────────────────────────────────────────────
    htl_path = os.path.join(output_dir, "epitope_report_htl.csv")
    htl_rows = [_htl_row(i, ep) for i, ep in enumerate(htl_epitopes, 1)]
    _write_csv(htl_path, HTL_COLUMNS, htl_rows)
    paths["htl_csv"] = htl_path

    # ── 3. B-cell CSV ─────────────────────────────────────────────────────────
    bcell_path = os.path.join(output_dir, "epitope_report_bcell.csv")
    bcell_rows = [_bcell_row(i, ep) for i, ep in enumerate(bcell_epitopes, 1)]
    _write_csv(bcell_path, BCELL_COLUMNS, bcell_rows)
    paths["bcell_csv"] = bcell_path

    # ── 4. Combined CSV (uses ALL_COLUMNS superset) ───────────────────────────
    all_path = os.path.join(output_dir, "epitope_report_all.csv")
    all_rows = []
    for row in ctl_rows:
        all_rows.append({col: row.get(col, "") for col in ALL_COLUMNS})
    for row in htl_rows:
        all_rows.append({col: row.get(col, "") for col in ALL_COLUMNS})
    for row in bcell_rows:
        all_rows.append({col: row.get(col, "") for col in ALL_COLUMNS})
    _write_csv(all_path, ALL_COLUMNS, all_rows)
    paths["all_csv"] = all_path

    # ── 5. JSON report (for frontend) ─────────────────────────────────────────
    json_path = os.path.join(output_dir, "epitope_report.json")
    report = {
        "generated_at": timestamp,
        "combined_ctl_population_coverage": round(combined_ctl_coverage * 100, 2),
        "summary": {
            "ctl_count": len(ctl_epitopes),
            "htl_count": len(htl_epitopes),
            "bcell_count": len(bcell_epitopes),
            "total": len(ctl_epitopes) + len(htl_epitopes) + len(bcell_epitopes),
        },
        "ctl_epitopes": ctl_rows,
        "htl_epitopes": htl_rows,
        "bcell_epitopes": bcell_rows,
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    paths["json"] = json_path

    print(f"  [OK]   Epitope report exported:")
    print(f"         CTL:    {ctl_path}  ({len(ctl_rows)} epitopes, 7 criteria)")
    print(f"         HTL:    {htl_path}  ({len(htl_rows)} epitopes, 10 criteria)")
    print(f"         B-cell: {bcell_path}  ({len(bcell_rows)} epitopes)")
    print(f"         All:    {all_path}")
    print(f"         JSON:   {json_path}")

    return {"paths": paths, "report": report}


def _write_csv(filepath: str, columns: list, rows: list):
    """Write rows to a CSV file."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
