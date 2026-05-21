"""structure_api.py
Lightweight RNA folding API for Module 2.
Prefers RNAfold (ViennaRNA). Falls back to a config-driven heuristic MFE estimate.
All magic numbers are read from config["module2"]["rnafold_heuristic"].
"""
import os
import re
import shutil
import subprocess

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

app = FastAPI()

# Load config for heuristic fallback coefficients
try:
    from default_config import DEFAULT_CONFIG as _CFG
    _HEURISTIC = _CFG.get("module2", {}).get("rnafold_heuristic", {})
except ImportError:
    _HEURISTIC = {}

# Heuristic coefficients — resolved from config with sensible defaults
_MFE_PER_NT       = _HEURISTIC.get("mfe_per_nt",        -0.25)
_GC_DEVIATION_MID = _HEURISTIC.get("gc_deviation_mid",   0.50)
_GC_WEIGHT        = _HEURISTIC.get("gc_weight",           1.5)


class FoldRequest(BaseModel):
    sequence: str


def _heuristic_mfe(seq: str) -> float:
    """
    Config-driven heuristic MFE estimate when RNAfold is unavailable.
    Formula: mfe_per_nt * len * (1 + gc_weight * (gc_ratio - gc_deviation_mid))
    Coefficients come from config["module2"]["rnafold_heuristic"].
    """
    gc_count = seq.count("G") + seq.count("C")
    gc_ratio = gc_count / len(seq)
    mfe = _MFE_PER_NT * len(seq) * (1.0 + _GC_WEIGHT * (gc_ratio - _GC_DEVIATION_MID))
    return round(mfe, 2)


@app.post("/api/fold")
def fold_sequence(body: FoldRequest):
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence")

    # ── Prefer real RNAfold ───────────────────────────────────────────────────
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
            if len(out) < 2:
                raise ValueError("Unexpected RNAfold output")

            m = re.match(r"([\.\(\)]+)\s+\(([-0-9\.]+)\)", out[1].strip())
            if not m:
                raise ValueError("Could not parse RNAfold output")

            return {
                "structure": m.group(1),
                "mfe":       float(m.group(2)),
                "method":    "rnafold",
            }

        except subprocess.CalledProcessError as e:
            raise HTTPException(
                status_code=500,
                detail=f"RNAfold failed: {e.stderr.decode('utf-8')}",
            )
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    # ── Heuristic fallback ────────────────────────────────────────────────────
    return {
        "structure": "." * len(seq),
        "mfe":       _heuristic_mfe(seq),
        "method":    "heuristic_estimate",
        "note": (
            f"RNAfold not found. MFE estimated via heuristic "
            f"(mfe_per_nt={_MFE_PER_NT}, gc_weight={_GC_WEIGHT}). "
            f"Install ViennaRNA for real predictions."
        ),
    }


if __name__ == "__main__":
    import uvicorn
    from fastapi.staticfiles import StaticFiles
    from fastapi.responses import FileResponse
    from fastapi.middleware.cors import CORSMiddleware
    
    # We should re-add the static files serving from the original file since user's snippet 
    # seems to have dropped it by only showing the API part mostly. Wait, looking at the code, 
    # the original file had CORS and static files logic, let's restore it to ensure it still serves the frontend.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
    )
    MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
    app.mount("/css", StaticFiles(directory=os.path.join(MODULE_DIR, "css")), name="css")
    app.mount("/js", StaticFiles(directory=os.path.join(MODULE_DIR, "js")), name="js")
    app.mount("/data", StaticFiles(directory=os.path.join(MODULE_DIR, "data")), name="data")

    @app.get("/")
    def serve_index():
        return FileResponse(os.path.join(MODULE_DIR, "index.html"))

    host = os.getenv("STRUCTURE_HOST", "0.0.0.0")
    try:
        port = int(os.getenv("STRUCTURE_PORT", 8000))
    except ValueError:
        port = 8000
    uvicorn.run(app, host=host, port=port)
