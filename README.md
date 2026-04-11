# 🧬 Research-Grade Vaccine Design & mRNA Optimization Pipeline

An end-to-end, research-grade platform for the computational design of multi-epitope vaccines and highly stable mRNA constructs. This pipeline integrates real-world biological constraints, population-scale demographics, and structural thermodynamics.

---

## 🚀 Research-Grade Highlights

- **Structural Selection Loop:** Optimization now generates multiple candidates and selects the one with the highest **$\Delta G$ stability**.
- **Population Analytics:** Research-grade heuristic mapping for global HLA demographics.
- **Unified Optimization:** GC-aware codon selection optimized for high-expression constructs.
- **Research Safety Floor:** Strict enforcement of a **0.75 Codon Adaptation Index (CAI)** target for elite constructs.
- **Adaptive Fallbacks:** Heuristic MFE estimation when `RNAfold` structural analysis is unavailable.
- **7-Step Bio-Safety Pipeline:** Detection and removal of cryptic splice sites, IRES motifs, and GC-rich hairpins.

---

## 📂 Project Architecture

```text
.
├── module1/                # Biological Analytics (Python)
│   ├── main.py             # Global orchestrator with population scaling
│   ├── default_config.py   # Centralized prediction & scoring weights
│   ├── advanced_filters.py # Toxicity & Allergenicity filters
│   └── scoring.py          # Multi-parameter affinity scoring
└── module2/                # mRNA Sequence Generator (JS/Python)
    ├── index.html          # Research Design dashboard
    ├── structure_api.py    # FastAPI backend for mRNA thermodynamics (with fallback)
    ├── js/config.js        # Centralized constants & data paths
    └── js/                 # Structural optimization engine
```

---

## ⚙️ Configuration & Environment

The platform is designed to be fully parameter-driven.

### Module 1: Python Parameters (`module1/default_config.py`)
| Key | Description | Default |
| :--- | :--- | :--- |
| `mhc1_method` | Prediction model for MHC-I | `netmhcpan` |
| `mhc2_method` | Prediction model for MHC-II | `netmhciipan` |
| `pop_coverage_bonus` | Weight boost for high-coverage epitopes | `5.0` |
| `linkers` | Linker strings for fusion | `AAY`, `GPGPG`, `KK` |

### Module 2: Environment Variables
- `STRUCTURE_HOST`: Interface to bind (default: `0.0.0.0`)
- `STRUCTURE_PORT`: Port to listen (default: `8000`)

---

## 🔬 Module 1: Biological Analytics
Module 1 identifies and selects high-affinity epitopes while enforcing strict research-grade biological constraints.

**Usage:**
```bash
python module1/main.py
```

---

## 💊 Module 2: Sequence Generator
Module 2 converts the antigen protein into a high-expression, stable mRNA sequence via structural selection logic.

### Research-Grade Loop:
1. **Multi-Candidate Generation:** Generates 3 unique optimized sequences per strategy.
2. **Structural Selection:** Computes $\Delta G$ for each and selects the most thermodynamically stable construct.
3. **Threshold Validation:** Verifies the final construct against `CAI > 0.75` and `GC 45-60%`.

**Setup Structural Backend:**
```bash
# Requires Python 3.x. Fallback MFE estimate active if RNAfold is missing.
python module2/structure_api.py
```

---

## 🛠 Technology Stack
- **Backend:** Python 3.x, FastAPI, Pandas.
- **Frontend:** Vanilla JavaScript (ES6+), Modern CSS.
- **Bio-Integration:** IEDB APIs, ToxinPred3 interface, ViennaRNA (RNAfold) bridge.

---

## 📜 License
*Research-Grade Development v4.0 - Optimized for Research and Vaccine Development.*