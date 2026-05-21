# 🧬 Research-Grade Vaccine Design & mRNA Optimization Pipeline

An end-to-end, research-grade platform for the computational design of multi-epitope vaccines and highly stable mRNA constructs. This pipeline integrates real-world biological constraints, structural thermodynamics, molecular docking, and cytokine predictions in a unified single-page application.

---

## 🚀 Research-Grade Highlights

- **Unified Dashboard:** Seamless one-click pipeline bridging protein antigen input to fully validated mRNA sequence, 3D structure, and docking analyses.
- **Structural Selection Loop:** Generates multiple mRNA candidates and selects the one with the highest **$\Delta G$ stability** using RNAfold or a deterministic MFE heuristic.
- **3D Structure Prediction:** Programmatic integration with Meta's **ESMFold API** for rapid tertiary structure modeling, complete with per-residue pLDDT confidence scoring and Ramachandran/Secondary structure insights.
- **Molecular Docking Simulation:** Auto-fetches standard receptors from RCSB (TLR4, TLR2, HLA-A*02:01) and estimates binding interactions with HDOCK submission capabilities.
- **Cytokine Induction Profiling:** Predicts Th1 (IFN-γ), Th2 (IL-4), and immunomodulatory (IL-10) responses using validated amino acid composition heuristics.
- **Real Population Analytics:** Mathematical probabilistic modeling of global demographics via IEDB exact logic.
- **Strict Hydrophobicity Control:** Heavy GRAVY scoring penalties enforcing structurally soluble constructs.
- **Cross-Strain Conservancy:** Sliding-window sequence alignment allowing biological motif variance.

---

## 📂 Project Architecture

```text
.
├── module1/                # Biological Analytics, Pipeline APIs & UI
│   ├── main.py             # Global orchestrator for epitope design
│   ├── module1_api.py      # Main FastAPI server (Unified Dashboard & Endpoints)
│   ├── index.html          # Unified single-page pipeline interface
│   ├── structure_3d.py     # ESMFold API integration for 3D PDB generation
│   ├── docking.py          # Receptor fetching and docking heuristics
│   ├── cytokine.py         # Th1/Th2 cytokine prediction tools
│   ├── default_config.py   # Centralized prediction & scoring weights
│   ├── advanced_filters.py # Toxicity & Allergenicity filters
│   └── scoring.py          # Multi-parameter affinity scoring
└── module2/                # mRNA Optimization & Folding Backend
    ├── structure_api.py    # FastAPI backend for RNAfold / MFE estimations
    ├── js/app.js           # Structural optimization logic module
    └── js/config.js        # Codon usage and threshold configurations
```

---

## ⚙️ How to Run the Unified Pipeline

The pipeline is fully automated and runs entirely locally across two inter-communicating API servers.

### Step 1: Start Module 2 (mRNA Structure Backend)
This server handles the structural evaluation ($\Delta G$) of generated mRNA sequences.
```bash
cd module2
python structure_api.py
```
*Runs on `http://localhost:8000`*

### Step 2: Start Module 1 (Pipeline API & Dashboard)
This server orchestrates the biological design, 3D folding, cytokines, and molecular docking.
```bash
cd module1
python module1_api.py
```
*Runs on `http://localhost:8080`*

### Step 3: Launch the Dashboard
Open your web browser and navigate to:
👉 **http://localhost:8080**

1. Enter your target Antigen sequence.
2. Adjust GC% and Stop Codon preferences if needed.
3. Click **"Run Full Pipeline"** and watch the automated workflow progress through all 6 phases.

---

## 🔬 Pipeline Phases

### Phase 1: Epitope Prediction (Module 1)
- Extracts MHC-I, MHC-II, and linear B-cell epitopes via IEDB APIs.
- Assesses antigenicity (VaxiJen logic), cross-strain conservancy, and toxicity.

### Phase 2: Vaccine Construct Assembly
- Filters epitopes based on affinity and population coverage.
- Joins epitopes using optimized linkers (`AAY` for CTL, `GPGPG` for HTL, `KK` for B-cell).
- Appends adjuvants (e.g., 50S ribosomal protein) and a 6xHis-tag for purification.

### Phase 3: Codon Optimization (Module 2)
- Reverse translates the protein sequence using species-specific codon usage tables (Human).
- Applies GC-balancing, removes forbidden motifs (e.g., poly-A, cryptic splice sites).
- Validates the Codon Adaptation Index (`CAI > 0.75`).

### Phase 4: 3D Structure Prediction
- Submits the vaccine construct to **ESMFold**.
- Returns a downloadable `.pdb` file and evaluates folding confidence via mean pLDDT scores.

### Phase 5: Cytokine Prediction
- Evaluates the construct's epitopes against established composition rules for **IFN-γ** (Th1), **IL-4** (Th2), and **IL-10** induction.

### Phase 6: Molecular Docking
- Downloads innate and adaptive immune receptors (TLR4, TLR2, HLA).
- Provides estimated binding energies and builds one-click submission links for HDOCK and PatchDock.

---

## 🛠 Technology Stack
- **Backend:** Python 3.x, FastAPI, Pandas, Requests.
- **Frontend:** HTML5, CSS3, Vanilla JavaScript (ES6+), HTML Canvas (RNA Secondary Structure Arc Diagrams).
- **Bio-Integration:** IEDB APIs, ESMFold API, RCSB PDB, ViennaRNA (RNAfold).

---

## 📜 License
*Research-Grade Development v4.1 - Optimized for Automated in silico Vaccine Design.*