from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import uvicorn
import os
import json
import requests as _requests

# Core pipeline modules
from main          import run_module1
from structure_3d  import predict_structure
from cytokine      import predict_cytokines, summarise_cytokines
from docking       import run_docking_analysis
from expression_vector import design_expression_vector
from secondary_structure import analyze_secondary_structure

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

class PipelineRequest(BaseModel):
    sequence: str

@app.post("/api/run")
def run_pipeline(body: PipelineRequest):
    """Module 1 only — returns the vaccine construct + epitope report."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence provided.")
    try:
        result = run_module1(seq)
        epitopes = result.get("top_epitopes", [])
        physchem = result.get("physicochemical", {})
        antigenicity = result.get("construct_antigenicity", {})
        return {
            "status": "success",
            "antigen_sequence": result.get("antigen_sequence", ""),
            "epitope_count": len(epitopes),
            "ctl_count": sum(1 for e in epitopes if e.get("type") == "MHC-I"),
            "htl_count": sum(1 for e in epitopes if e.get("type") == "MHC-II"),
            "bcell_count": sum(1 for e in epitopes if e.get("type") == "B-cell"),
            "molecular_weight_kda": physchem.get("molecular_weight_kda"),
            "antigenicity_score": antigenicity.get("score"),
            "is_antigen": antigenicity.get("is_antigen"),
            "combined_ctl_coverage": result.get("combined_ctl_coverage", 0),
            "epitope_report": result.get("epitope_report", {}),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/epitope-report")
def get_epitope_report():
    """Return the latest epitope selection report as JSON."""
    json_path = os.path.join(MODULE_DIR, "output", "epitope_report.json")
    if not os.path.exists(json_path):
        raise HTTPException(status_code=404,
                            detail="No epitope report found. Run the pipeline first.")
    with open(json_path, encoding="utf-8") as f:
        report = json.load(f)
    return report


@app.get("/api/epitope-csv/{report_type}")
def download_epitope_csv(report_type: str):
    """
    Download epitope report as CSV.
    report_type: 'ctl' | 'htl' | 'bcell' | 'all'
    """
    valid = {"ctl", "htl", "bcell", "all"}
    if report_type not in valid:
        raise HTTPException(status_code=400,
                            detail=f"Invalid type '{report_type}'. Use: {', '.join(valid)}")

    csv_path = os.path.join(MODULE_DIR, "output", f"epitope_report_{report_type}.csv")
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404,
                            detail=f"No {report_type} CSV found. Run the pipeline first.")

    with open(csv_path, encoding="utf-8") as f:
        content = f.read()

    filename = f"epitope_report_{report_type}.csv"
    return Response(
        content=content,
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


@app.post("/api/run-full")
def run_full_pipeline(body: PipelineRequest):
    """Full pipeline: Module 1 → Module 2. Returns vaccine construct + mRNA metrics."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence provided.")
    
    try:
        result = run_module1(seq)
        antigen = result.get("antigen_sequence", "")
        if not antigen:
            raise ValueError("Module 1 produced no antigen sequence.")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Module 1 failed: {e}")

    rna_seq = antigen.replace("T", "U")
    mod2_result = {"mfe": None, "method": "not_available", "structure": None}

    try:
        import sys, os
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "module2"))
        from structure_api import fold_sequence, FoldRequest
        res = fold_sequence(FoldRequest(sequence=rna_seq))
        mod2_result = res if isinstance(res, dict) else res.dict()
        mod2_result["port_used"] = "internal"
    except Exception as ex:
        print("Internal fold failed, trying HTTP fallback:", ex)
        for port in [8001, 8000]:  # Only attempt separate background server ports
            try:
                resp = _requests.post(
                    f"http://localhost:{port}/api/fold",
                    json={"sequence": rna_seq},
                    timeout=5
                )
                if resp.status_code == 200:
                    mod2_result = resp.json()
                    mod2_result["port_used"] = port
                    break
            except Exception:
                continue


    return {
        "status": "success",
        "antigen_sequence": antigen,
        "combined_ctl_coverage": result.get("combined_ctl_coverage", 0),
        "epitope_report": result.get("epitope_report", {}),
        "mrna": {
            "rna_sequence": rna_seq[:100] + "..." if len(rna_seq) > 100 else rna_seq,
            "mfe_kcal_mol": mod2_result.get("mfe"),
            "structure_method": mod2_result.get("method"),
            "module2_available": mod2_result.get("mfe") is not None,
        }
    }


class EpitopeListRequest(BaseModel):
    epitopes: list  # list of epitope dicts from Module 1


@app.post("/api/structure")
def get_structure(body: PipelineRequest):
    """Priority 1: ESMFold 3D structure prediction."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence.")
    try:
        result = predict_structure(seq, output_dir="output")
        result.pop("pdb_content", None)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/pdb")
def download_pdb():
    """Download the latest generated PDB file."""
    pdb_path = os.path.join("output", "structure.pdb")
    if not os.path.exists(pdb_path):
        raise HTTPException(status_code=404, detail="No PDB file found. Run structure prediction first.")
    with open(pdb_path) as f:
        content = f.read()
    return Response(content=content, media_type="chemical/x-pdb",
                    headers={"Content-Disposition": "attachment; filename=vaccine_structure.pdb"})


@app.post("/api/cytokines")
def get_cytokines(body: EpitopeListRequest):
    """Priority 2: IFN-γ, IL-4, IL-10 prediction for all epitopes."""
    epitopes = body.epitopes
    if not epitopes:
        raise HTTPException(status_code=400, detail="No epitopes provided.")
    try:
        enriched = predict_cytokines(epitopes)
        summary  = summarise_cytokines(enriched)
        return {"status": "success", "epitopes": enriched, "summary": summary}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/docking")
def get_docking(body: PipelineRequest):
    """Priority 3: Molecular docking vs TLR4, TLR2, HLA-A*02:01."""
    pdb_path = os.path.join("output", "structure.pdb")
    try:
        result = run_docking_analysis(pdb_path, output_dir="output")
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/vector")
def get_expression_vector(body: PipelineRequest):
    """§2.33-2.34: Expression vector design + in silico cloning."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence.")
    try:
        result = design_expression_vector(seq)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/secondary-structure")
def get_secondary_structure(body: PipelineRequest):
    """SOPMA & PSIPRED secondary structure prediction and report generation."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence.")
    try:
        result = analyze_secondary_structure(seq, output_dir=OUTPUT_DIR)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/secondary-structure-csv")
def download_secondary_structure_csv():
    """Download per-residue SOPMA & PSIPRED secondary structure CSV report."""
    csv_path = os.path.join(OUTPUT_DIR, "secondary_structure_report.csv")
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail="No secondary structure CSV found. Run secondary structure analysis first.")
    with open(csv_path, encoding="utf-8") as f:
        content = f.read()
    return Response(content=content, media_type="text/csv",
                    headers={"Content-Disposition": "attachment; filename=secondary_structure_report.csv"})



MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
if os.environ.get("VERCEL"):
    OUTPUT_DIR = "/tmp/output"
else:
    OUTPUT_DIR = os.path.join(MODULE_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
app.mount("/output", StaticFiles(directory=OUTPUT_DIR), name="output")

@app.get("/")
def serve_index():
    return FileResponse(os.path.join(MODULE_DIR, "index.html"))

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080)
