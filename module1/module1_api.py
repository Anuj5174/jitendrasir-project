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
    """Module 1 only — returns the vaccine construct."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence provided.")
    try:
        result = run_module1(seq)
        return {
            "status": "success",
            "antigen_sequence": result.get("antigen_sequence", "")
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/run-full")
def run_full_pipeline(body: PipelineRequest):
    """Full pipeline: Module 1 → Module 2. Returns vaccine construct + mRNA metrics."""
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence provided.")
    
    # ── Step 1: Run Module 1 ──────────────────────────────────────────────────
    try:
        result = run_module1(seq)
        antigen = result.get("antigen_sequence", "")
        if not antigen:
            raise ValueError("Module 1 produced no antigen sequence.")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Module 1 failed: {e}")

    # ── Step 2: Trigger Module 2 fold API ─────────────────────────────────────
    rna_seq = antigen.replace("T", "U")  # DNA → RNA
    mod2_result = {"mfe": None, "method": "not_available", "structure": None}

    # Try module2 structure API (port 8001 or 8000)
    for port in [8001, 8000]:
        try:
            resp = _requests.post(
                f"http://localhost:{port}/api/fold",
                json={"sequence": rna_seq},
                timeout=30
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
        # Don't serialise raw PDB content in JSON — use /api/pdb endpoint
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


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

@app.get("/")
def serve_index():
    return FileResponse(os.path.join(MODULE_DIR, "index.html"))

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080)
