# 🧬 Enterprise Vaccine Design & mRNA Optimization Pipeline

An end-to-end, clinical-tier platform for the computational design of multi-epitope vaccines and highly stable mRNA constructs. This pipeline integrates real-world biological constraints, population-scale demographics, and unified structural thermodynamics.

---

## 🚀 Key Highlights

- **Clinical-Tier Screening:** Fixed integration with **NetMHCpan**, **NetMHCIIpan**, and **BepiPred** via IEDB APIs.
- **Unified Optimization:** GC-aware codon selection that ensures high expression without sacrificing **$\Delta G$ stability**.
- **Therapeutic Safety Floor:** Automated enforcement of a **0.7 Codon Adaptation Index (CAI)** threshold.
- **Population-Scale Analytics:** Integrated Global HLA Population Coverage and Conservancy metrics.
- **Context-Aware CPB:** Professional-grade Codon Pair Bias optimization using a curated human weight matrix.
- **7-Step Bio-Safety Pipeline:** Detection and removal of cryptic splice sites, IRES motifs, and GC-rich hairpins.

---

## 📂 Project Architecture

```text
.
├── module1/                # Clinical Epitope Selection (Python)
│   ├── main.py             # Global orchestrator with population scaling
│   ├── default_config.py   # Centralized prediction & scoring weights
│   ├── advanced_filters.py # Toxicity (ToxinPred3) & Allergenicity filters
│   └── scoring.py          # Multi-parameter affinity scoring
└── module2/                # mRNA Sequence Optimizer (JS/Python)
    ├── index.html          # Therapeutic Design dashboard
    ├── structure_api.py    # FastAPI backend for mRNA thermodynamics
    ├── js/config.js        # Centralized web UI constants & data paths
    └── js/                 # Unified GC-Codon optimization engine
```

---

## ⚙️ Configuration & Environment

The platform is designed to be fully parameter-driven without hardcoded values.

### Module 1: Python Parameters (`module1/default_config.py`)
| Key | Description | Default |
| :--- | :--- | :--- |
| `mhc1_method` | Prediction model for MHC-I | `netmhcpan` |
| `mhc2_method` | Prediction model for MHC-II | `netmhciipan` |
| `pop_coverage_bonus` | Weight boost for high-coverage epitopes | `5.0` |
| `linkers` | Linker strings for fusion (MHC-I/II/B) | `AAY`, `GPGPG`, `KK` |

### Module 2: Environment Variables
The structural backend supports standard environment variables for deployment:
- `STRUCTURE_HOST`: Interface to bind (default: `0.0.0.0`)
- `STRUCTURE_PORT`: Port to listen (default: `8000`)

---

## 🔬 Module 1: Epitope Selection (Clinical Tier)
Module 1 identifies and selects high-affinity epitopes while enforcing strict biological and demographic constraints.

**Usage:**
```bash
python module1/main.py
```
*Input secondary protein sequences via CLI or pipe.*

---

## 💊 Module 2: mRNA Optimizer (Therapeutic Tier)
Module 2 converts the antigen protein into a high-expression, stable mRNA sequence via a unified GC-Codon optimization strategy.

### The Unified Strategy:
1. **GC-Aware Translation:** Codons are selected based on $Score = Frequency \times CPB \times GC\_Bias$.
2. **Safe GC Refinement:** Synonymous swaps are performed only if local windows deviate from 52%, provided CAI stays above 0.7.
3. **Folding Thermodynamics:** Real-time calculation of $ \Delta G $ via the FastAPI backend.

**Setup Structural Backend:**
```bash
# Requires Python 3.x and ViennaRNA (RNAfold) in system path
python module2/structure_api.py
```

**Launch Dashboard:**
Serve the root directory via any local web server (e.g., `python -m http.server 3000`) and open `http://localhost:3000/module2/`.

---

## 🛠 Technology Stack
- **Backend:** Python 3.x, FastAPI, Pandas, Subprocess logic.
- **Frontend:** Vanilla JavaScript (ES6+), Modern CSS3 (Glassmorphism).
- **Bio-Integration:** IEDB APIs, ToxinPred3 interface, ViennaRNA (RNAfold) structural bridge.

---

## 📜 License
*Enterprise Clinical Tier v3.0 - Optimized for Clinical Research and Vaccine Development.*