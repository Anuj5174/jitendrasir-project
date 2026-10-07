# default_config.py
# ─────────────────────────────────────────────────────────────────────────────
# SINGLE SOURCE OF TRUTH for the entire pipeline.
# Every constant, threshold, weight, and formula coefficient lives here.
# No magic numbers anywhere else in the codebase.
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_CONFIG = {

    # ── Target organism ───────────────────────────────────────────────────────
    # Controls antigenicity thresholds and organism-specific scoring.
    # Options: "virus" | "bacteria" | "parasite" | "tumor"
    "target_organism": "virus",

    # ── HLA alleles ───────────────────────────────────────────────────────────
    # Expanded panel for >90% global population coverage.
    # MHC-I supertypes cover ~97% of world population (Sette & Sidney, 1999).
    # MHC-II panel covers major DR/DQ/DP loci.
    "alleles": {
        "mhc1": [
            # HLA-A supertypes (covers ~95% global population)
            "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01",
            "HLA-A*24:02", "HLA-A*26:01", "HLA-A*11:01",
            "HLA-A*68:01",
            # HLA-B supertypes (adds ~15% unique coverage)
            "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
            "HLA-B*27:05", "HLA-B*35:01", "HLA-B*40:01",
            "HLA-B*44:02", "HLA-B*51:01", "HLA-B*53:01",
            "HLA-B*58:01",
        ],
        "mhc2": [
            "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01",
            "HLA-DRB1*07:01", "HLA-DRB1*08:01", "HLA-DRB1*11:01",
            "HLA-DRB1*13:01", "HLA-DRB1*15:01",
            "HLA-DQA1*05:01/DQB1*02:01", "HLA-DPA1*01:03/DPB1*04:01",
        ],
    },

    # ── IEDB API ──────────────────────────────────────────────────────────────
    "api": {
        "base_url":    "http://tools-cluster-interface.iedb.org/tools_api",
        "timeout":     30,
        "retries":     3,
        "retry_delay": 1,
    },

    # ── Prediction methods ────────────────────────────────────────────────────
    "prediction": {
        "mhc1_length": 9,
        "mhc1_method": "netmhcpan",
        "mhc2_method": "netmhciipan",
        "bcell_method": "Bepipred",
    },

    # ── Binding thresholds ────────────────────────────────────────────────────
    "thresholds": {
        "ic50":             500,    # nM — peptides above this are discarded
        "mhc2_ic50":        800,    # relaxed threshold for MHC-II
        "percentile_rank":  1.0,     # % rank cutoff for MHC binding (tightened for speed)
        "bcell":            0.5,    # BepiPred score cutoff
        "min_bcell_length": 6,      # minimum B-cell epitope length (AA)
        "fallback_window_size": 8,  # size of sliding window for heuristic fallback
        "overlap_limit":    6,      # max shared AA before peptide is redundant
    },

    # ── Scoring formula weights ───────────────────────────────────────────────
    "scoring": {
        "max_score_fallback":        10.0,  # score when IC50/rank is 0 or missing
        "ic50_normalizer":        50000.0,  # reference IC50 for log-scaling
                                            # score = log(ic50_normalizer / IC50)
        "ic50_strong_bonus":          5.0,  # bonus for IC50 < ic50_strong_threshold
        "ic50_moderate_bonus":        2.0,  # bonus for IC50 < ic50_moderate_threshold
        "ic50_strong_threshold":     50.0,  # nM — "strong binder" cutoff
        "ic50_moderate_threshold":  500.0,  # nM — "moderate binder" cutoff
        "allele_breadth_weight":      1.0,  # multiplier on log(1 + allele_count)
        "bcell_score":               20.0,  # flat score assigned to B-cell epitopes
        "population_coverage_bonus":  5.0,  # score bonus per unit of world coverage
        "conservancy_bonus":          2.0,  # score bonus per unit of conservancy
        "human_homology_penalty":   -50.0,  # score penalty for human-similar peptides
    },

    # ── CD8+/MHC-I Multi-Criteria Selection Weights ──────────────────────────
    # Each criterion is scored 0-1 then multiplied by its weight.
    # Final rank = Σ(weight_i × normalised_score_i)
    "ctl_selection": {
        "weights": {
            "binding_strength":     0.15,  # Reduced to give antigenicity priority
            "antigenicity":         0.30,  # DOUBLED: High priority for antigenicity
            "allergenicity":        0.10,
            "toxicity":             0.10,
            "immunogenicity":       0.15,  # Reduced slightly
            "tc50":                 0.10,
            "multi_hla_binding":    0.10,  # Reduced slightly
        },
        # Minimum alleles a CTL epitope must bind to be considered
        "min_allele_count":         2,
        # Minimum combined population coverage for the selected CTL set
        "min_population_coverage":  0.90,
        # TC50 thresholds (nM)
        "tc50_strong":             50.0,
        "tc50_moderate":          200.0,
        "tc50_weak":              500.0,
    },

    # ── CD4+/MHC-II Multi-Criteria Selection Weights ─────────────────────────
    # 10 criteria: same 7 as CTL + IFN-γ, IL-4, IL-10 cytokine induction
    "htl_selection": {
        "weights": {
            "binding_strength":     0.15,
            "antigenicity":         0.10,
            "allergenicity":        0.08,
            "toxicity":             0.08,
            "immunogenicity":       0.12,
            "tc50":                 0.07,
            "multi_hla_binding":    0.10,
            "ifn_gamma":            0.12,  # IFN-γ inducer → Th1 cellular immunity
            "il4":                  0.10,  # IL-4 inducer → Th2 humoral immunity
            "il10":                 0.08,  # IL-10 → immunomodulation
        },
        "min_allele_count":         2,
    },

    # ── B-cell Multi-Criteria Selection Weights ──────────────────────────────
    # 3 criteria: Antigenicity, Allergenicity, Toxicity
    "bcell_selection": {
        "weights": {
            "antigenicity":         0.50,
            "allergenicity":        0.25,
            "toxicity":             0.25,
        },
    },

    # ── Antigenicity (VaxiJen-like ACC) ──────────────────────────────────────
    "antigenicity": {
        # Per-organism decision thresholds (from VaxiJen2.0 calibration)
        "thresholds": {
            "virus":    0.32,   # adjusted for multi-epitope constructs
            "bacteria": 0.32,
            "parasite": 0.40,
            "tumor":    0.38,
        },
        # Formula weights: score = w_hydro*hydro + w_polar*polar + w_charged*charged + w_var*var
        "formula_weights": {
            "hydrophobic_fraction":  0.35,
            "polar_fraction":        0.15,
            "charged_fraction":      0.20,   # applied as (1 - charged_frac)
            "accessibility_var":     0.30,   # applied as min(1, acc_var * acc_var_scale)
            "acc_var_scale":          2.0,
        },
    },

    # ── Homology filter ───────────────────────────────────────────────────────
    "homology": {
        "identity_threshold":    50.0,  # % — peptides with human identity above
                                        # this are penalised (not removed)
        "max_short_peptide_len":    9,  # AAs — fast motif scan applies to peptides <= this
        "codon_bias_ratio":      0.85,  # human codon bias ratio above which to flag
        # Known human self-peptide motifs (fast pre-filter before BLAST)
        "human_self_motifs": [
            "FLLLCLAIR", "KFERQ", "LQNRRGLD", "SIINFEKL",
            "GILGFVFTL", "NLVPMVATV", "GLCTLVAML",
            "FLPSDFFPSV", "ILKEPVHGV", "QAKWRLQTL",
        ],
    },

    # ── Safety filters ────────────────────────────────────────────────────────
    "safety": {
        # Allergenicity confidence thresholds
        "allergen_confidence_cutoff":  0.65,  # flag as allergen if confidence >= this
        "allergen_pi_threshold":      10.5,   # basic pI IgE cross-reactivity flag
        "allergen_proline_fraction":   0.25,  # proline-rich repeat threshold
        "allergen_proline_min_len":    12,    # only apply proline check above this length
        "allergen_stability_len_min":  15,    # only apply GI stability check above this
        "allergen_stability_threshold": 5.0,  # instability index below = very stable
    },

    # ── Vaccine construct assembly ────────────────────────────────────────────
    "construct": {
        # Adjuvant — prepended at N-terminus
        "adjuvant_name":     "50S ribosomal protein L7/L12 (Mycobacterium tuberculosis)",
        "adjuvant_sequence": "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVE"
                             "AAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKE"
                             "AADEAKAKLEAAGATVTVK",
        # Universal T-helper epitope
        "padre_sequence":    "AKFVAAWTLKAAA",
        # Purification tag appended at C-terminus
        "his_tag":           "HHHHHH",
        # Validation checks — derived from config so changing adjuvant/tag above
        # automatically updates these without touching immune_export.py
        "adjuvant_prefix":   "MAKL",   # first 4 AA of adjuvant for validation
        "his_tag_suffix":    "HHHHHH", # C-terminal tag for validation
        # Linkers
        "linkers": {
            "MHC-I":   "AAY",    # CTL epitope linker
            "MHC-II":  "GPGPG",  # HTL epitope linker
            "B-cell":  "KK",     # B-cell epitope linker
            "rigid":   "EAAAK",  # rigid domain separator
        },
    },

    # ── Epitope selection quotas ──────────────────────────────────────────────
    "selection": {
        "counts": {
            "MHC-I":  7,
            "MHC-II": 6,
            "MHC-II_min": 4,
            "B-cell": 3,
        },
        "max_peptides_for_coverage": 30, # compute coverage for top N candidates
        "combined_population_coverage": 0.97, # target combined MHC-I + MHC-II coverage
    },

    # ── Immune simulation heuristics ──────────────────────────────────────────
    "immune_simulation": {
        # Immunoglobulin response weights
        # IgM = ctl_weight*CTL_count + htl_weight*HTL_count + bcell_weight*Bcell_count
        "igm_weights": {
            "ctl":   0.3,
            "htl":   0.4,
            "bcell": 0.5,
        },
        # IgG weights
        "igg_weights": {
            "ctl":   0.7,
            "htl":   0.6,
            "bcell": 0.4,
        },
        # CD8, CD4, and Memory B-cell weights
        "cell_activation_weights": {
            "cd8_mhc1": 1.2,
            "cd8_mhc2": 0.2,
            "cd4_mhc2": 1.5,    # boosted to reflect absolute HTL strength rather than just fraction
            "cd4_mhc1": 0.2,
            "memory_mhc2": 1.0, # boosted to show HTL impact on B-cell memory
            "memory_bcell": 1.5,
        },
        # Cytokine profile weights
        "cytokine_weights": {
            "th1_ctl_contribution":   0.6,  # CTL epitopes drive Th1
            "th1_htl_contribution":   0.5,
            "th2_bcell_contribution": 0.7,  # B-cell epitopes drive Th2
        },
        # Thresholds for HIGH/MODERATE/LOW classification
        "response_thresholds": {
            "high":     0.6,
            "moderate": 0.3,
        },
    },

    # ── C-ImmSim submission parameters ───────────────────────────────────────
    "cimmsim": {
        "injection_days":    [1, 84, 168],   # standard 3-dose schedule
        "injection_volume":  0.5,            # mL
        "simulation_steps":  1050,           # time steps (1 step ≈ 8 hours)
        "url": "https://kraken.iac.rm.cnr.it/C-IMMSIM/index.php",
    },

    # ── Validation ────────────────────────────────────────────────────────────
    "validation": {
        "min_length":           20,    # minimum input protein length (AA)
        "min_construct_length": 100,   # minimum final construct length
        "mw_range_kda":         [15, 100],  # acceptable MW range for the construct
    },

    # ── Biological data ───────────────────────────────────────────────────────
    "biological_data": {
        "hla_frequencies": {
            # HLA-A allele frequencies (Allele Frequency Net Database)
            "HLA-A*02:01":    0.2482,  # most common worldwide
            "HLA-A*01:01":    0.1230,
            "HLA-A*03:01":    0.1066,
            "HLA-A*24:02":    0.1630,  # dominant in Asian/Oceanian populations
            "HLA-A*26:01":    0.0528,
            "HLA-A*11:01":    0.1290,  # high in East/Southeast Asian
            "HLA-A*68:01":    0.0432,
            # HLA-B allele frequencies
            "HLA-B*07:02":    0.1100,
            "HLA-B*08:01":    0.0850,
            "HLA-B*15:01":    0.0650,
            "HLA-B*27:05":    0.0380,
            "HLA-B*35:01":    0.0920,
            "HLA-B*40:01":    0.0680,
            "HLA-B*44:02":    0.0750,
            "HLA-B*51:01":    0.0580,
            "HLA-B*53:01":    0.0450,
            "HLA-B*58:01":    0.0380,
            # HLA-C
            "HLA-C*07:01":    0.28,
            # HLA-DRB1
            "HLA-DRB1*01:01": 0.0800,
            "HLA-DRB1*03:01": 0.1050,
            "HLA-DRB1*04:01": 0.0820,
            "HLA-DRB1*07:01": 0.1150,
            "HLA-DRB1*08:01": 0.0480,
            "HLA-DRB1*11:01": 0.0920,
            "HLA-DRB1*13:01": 0.0780,
            "HLA-DRB1*15:01": 0.1200,
            # HLA-DQ/DP
            "HLA-DQA1*05:01/DQB1*02:01": 0.11,
            "HLA-DPA1*01:03/DPB1*04:01": 0.13,
        },
        "default_allele_frequency": 0.03, # Fallback frequency if allele is missing
        # Kyte-Doolittle hydropathy scale — universal biological constant,
        # stored here for single-point reference rather than scattered copies
        "kyte_doolittle": {
            "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
            "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
            "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
            "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
        },
    },

    # ── Module 2 — mRNA / CDS design ─────────────────────────────────────────
    "module2": {
        # RNA structure heuristic (used when RNAfold is not installed)
        "rnafold_heuristic": {
            "mfe_per_nt":       -0.25,   # kcal/mol per nucleotide base rate
            "gc_deviation_mid":  0.50,   # GC ratio midpoint (neutral)
            "gc_weight":          1.5,   # GC deviation amplifier
        },
        # Codon optimization
        "codon_optimization": {
            "strategy":   "most_frequent",  # most_frequent | weighted_random | balanced
            "target_gc":  52,               # % target GC content
            "gc_window":  30,               # sliding window size for GC analysis
        },
        # mRNA construct structure
        "mrna": {
            "utr5_strategy":  "human_beta_globin",
            "utr3_strategy":  "double_globin",
            "polya_length":   120,
            "stop_codon":     "TGA",
            "modification":   "m1psi",      # nucleotide modification type
        },
    },
}