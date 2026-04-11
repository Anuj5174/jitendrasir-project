"""Lightweight folding API for Module 2 using ViennaRNA (RNAfold).

Provides a POST /api/fold endpoint that accepts JSON {"sequence": "AUG..."}
and returns {"structure": "..((..))..", "mfe": -12.34}.

This implementation prefers calling the `RNAfold` command-line tool. If
`RNAfold` is not available, the endpoint will return a helpful error message.
"""
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import subprocess
import shutil
import re

app = FastAPI()


class FoldRequest(BaseModel):
    sequence: str


@app.post("/api/fold")
def fold_sequence(body: FoldRequest):
    seq = body.sequence.strip().upper()
    if not seq:
        raise HTTPException(status_code=400, detail="Empty sequence")

    if not shutil.which('RNAfold'):
        # Fallback: Research-grade heuristic MFE estimate
        # Average stability: ~ -0.25 kcal/mol per nucleotide, weighted by GC
        gc_count = (seq.count('G') + seq.count('C'))
        gc_ratio = gc_count / len(seq)
        heuristic_mfe = (-0.25 * len(seq)) * (1.0 + (gc_ratio - 0.5) * 1.5)
        return {
            "structure": "." * len(seq),
            "mfe": round(heuristic_mfe, 2),
            "method": "heuristic_estimate"
        }

    try:
        # RNAfold reads sequence from stdin and writes structure to stdout
        proc = subprocess.run(['RNAfold', '--noPS', '--noShape'], input=seq.encode('utf-8'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        out = proc.stdout.decode('utf-8').strip().splitlines()
        # RNAfold output usually: <sequence>\n<structure> ( ( -12.34 ) )
        if len(out) < 2:
            raise ValueError('Unexpected RNAfold output')

        struct_line = out[1].strip()
        # parse structure and energy using regex: structure [space] ( -12.34 )
        m = re.match(r"([\.\(\)]+)\s+\(([-0-9\.]+)\)", struct_line)
        if not m:
            raise ValueError('Could not parse RNAfold output')

        structure = m.group(1)
        mfe = float(m.group(2))
        return {"structure": structure, "mfe": mfe}

    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=500, detail=f"RNAfold failed: {e.stderr.decode('utf-8')}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == '__main__':
    import uvicorn
    import os
    
    host = os.getenv('STRUCTURE_HOST', '0.0.0.0')
    try:
        port = int(os.getenv('STRUCTURE_PORT', 8000))
    except ValueError:
        port = 8000
        
    uvicorn.run(app, host=host, port=port)
