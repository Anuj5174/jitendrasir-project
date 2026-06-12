"""structure_api.py
Lightweight RNA folding API + static file server for Module 2.
Prefers RNAfold (ViennaRNA). Falls back to a config-driven heuristic MFE estimate.
"""
import os
import re
import shutil
import subprocess

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

app = FastAPI()

# CORS — must be at app level so cross-origin JS imports from port 8080 work
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Load config for heuristic fallback coefficients
try:
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "module1"))
    from default_config import DEFAULT_CONFIG as _CFG
    _HEURISTIC = _CFG.get("module2", {}).get("rnafold_heuristic", {})
except Exception:
    _HEURISTIC = {}

# Heuristic coefficients — resolved from config with sensible defaults
_MFE_PER_NT       = _HEURISTIC.get("mfe_per_nt",        -0.25)
_GC_DEVIATION_MID = _HEURISTIC.get("gc_deviation_mid",   0.50)
_GC_WEIGHT        = _HEURISTIC.get("gc_weight",           1.5)


class FoldRequest(BaseModel):
    sequence: str


def _heuristic_mfe(seq: str) -> float:
    gc_count = seq.count("G") + seq.count("C")
    gc_ratio = gc_count / len(seq)
    mfe = _MFE_PER_NT * len(seq) * (1.0 + _GC_WEIGHT * (gc_ratio - _GC_DEVIATION_MID))
    return round(mfe, 2)


@app.post("/api/fold")
def fold_sequence(body: FoldRequest):
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence")

    # ── Method 1: ViennaRNA Python binding (pip install ViennaRNA) ────────────
    try:
        import RNA
        fc = RNA.fold_compound(seq)
        structure, mfe = fc.mfe()
        return {
            "structure": structure,
            "mfe":       round(mfe, 2),
            "method":    "rnafold_python",
            "note":      "Computed via ViennaRNA Python binding (pip)."
        }
    except ImportError:
        pass  # Fall through to next method

    # ── Method 2: RNAfold binary (if installed separately) ───────────────────
    if shutil.which("RNAfold"):
        try:
            proc = subprocess.run(
                ["RNAfold", "--noPS", "--noShape"],
                input=seq.encode("utf-8"),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            out = proc.stdout.decode("utf-8").strip().splitlines()
            if len(out) >= 2:
                m = re.match(r"([\.\(\)]+)\s+\(([-0-9\.]+)\)", out[1].strip())
                if m:
                    return {"structure": m.group(1), "mfe": float(m.group(2)), "method": "rnafold_binary"}
        except Exception:
            pass

    # ── Method 3: Heuristic fallback ─────────────────────────────────────────
    return {
        "structure": "." * len(seq),
        "mfe":       _heuristic_mfe(seq),
        "method":    "heuristic_estimate",
        "note":      "ViennaRNA not available. Install with: pip install ViennaRNA",
    }


# ── Serve Module 2 static files (JS, CSS, data) ───────────────────────────────
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

app.mount("/css",  StaticFiles(directory=os.path.join(MODULE_DIR, "css")),  name="css")
app.mount("/js",   StaticFiles(directory=os.path.join(MODULE_DIR, "js")),   name="js")
app.mount("/data", StaticFiles(directory=os.path.join(MODULE_DIR, "data")), name="data")


@app.get("/")
def serve_index():
    return FileResponse(os.path.join(MODULE_DIR, "index.html"))


if __name__ == "__main__":
    import uvicorn
    host = os.getenv("STRUCTURE_HOST", "0.0.0.0")
    try:
        port = int(os.getenv("STRUCTURE_PORT", 8000))
    except ValueError:
        port = 8000
    uvicorn.run(app, host=host, port=port)
