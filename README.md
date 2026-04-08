# Epitope Vaccine Construction Pipeline

An automated bioinformatics pipeline that transforms a viral or pathogenic protein sequence into a single, optimized multi-epitope vaccine construct ready for synthesis.

## 🎯 Objective
This system intelligently queries immune databases (IEDB), identifies the most potent immunogenic regions (MHC-I, MHC-II, and B-cell), mathematically scores them for efficacy, and safely fuses them into a single, highly efficient artificial antigen.

## ⚙️ Core Operations (The Pipeline)

1. **Intelligent Targeting**: Simultaneously interfaces with immune databases via the IEDB API to predict the strongest binding regions for **MHC-I** (CD8+ T-cells), **MHC-II** (CD4+ Helper T-cells), and **B-cells** (Antibodies).
2. **Advanced Scoring Algorithm**: Evaluates candidate peptides using a biological formula (`-log(IC50) + Allele Coverage`) to prioritize regions that are not only potent but also effective across a broad human population. 
3. **Diversity & Deduplication**: Employs a strict algorithmic filter to enforce structural diversity. Any candidate sharing a biological motif (6+ amino acid overlap) with an already-selected region is instantly rejected, ensuring the final vaccine spans distinct structural regions.
4. **Precision Assembly**: Extracts exactly **7 MHC-I**, **5 MHC-II**, and **3 B-cell** epitopes, seamlessly merging their junction boundaries (`merge_junction`) with structured biochemical linkers (`AAY`, `GPGPG`, `KK`) to avoid biological bloat or stacking.

## 🚀 Usage

```bash
# Run the pipeline
python main.py
```

Provide the tool with a raw protein sequence when prompted. The script will perform prediction, filtering, scoring, and fusion, outputting the final engineered, multi-faceted antigen sequence.

## 🛠️ Configuration
The entire pipeline is dynamically governed by `default_config.py`. Centralized parameters include:
- **Prediction Methods & Alleles:** Choose tools like BepiPred, NetMHCpan, NetMHCIIpan.
- **Scientific Thresholds:** Adjust `ic50`, percentile ranks, and the biological `overlap_limit`.
- **Selection Quotas:** Define exact target counts for MHC-I, MHC-II, and B-cell epitopes.
- **Linkers & Fallbacks:** Set biochemical linkers and algorithm scoring fallbacks.

## ✨ The Result
A highly engineered, non-redundant, and structurally clean antigen sequence perfectly balanced for comprehensive immune stimulation—all driven by a centralized configuration file.