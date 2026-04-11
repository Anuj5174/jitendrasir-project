DEFAULT_CONFIG = {
    "alleles": {
        "mhc1": ["HLA-A*02:01"],
        "mhc2": ["HLA-DRB1*01:01"]
    },

    "thresholds": {
        "ic50": 500,
        "percentile_rank": 10,
        "bcell": 0.5,
        "min_bcell_length": 6,
        "overlap_limit": 6
    },

    "prediction": {
        "mhc1_length": 9,
        "mhc1_method": "netmhcpan",
        "mhc2_method": "netmhciipan",
        "bcell_method": "bepipred"
    },

    "scoring": {
        "allele_weight": 2,
        "bcell_score": 20,
        "max_score_fallback": 10.0,
        "population_coverage_bonus": 5.0,
        "conservancy_bonus": 2.0
    },

    "fusion": {
        "linkers": {
            "MHC-I": "AAY",
            "MHC-II": "GPGPG",
            "B-cell": "KK"
        }
    },

    "selection": {
        "counts": {
            "MHC-I": 7,
            "MHC-II": 5,
            "B-cell": 3
        }
    },

    "api": {
        "base_url": "http://tools-cluster-interface.iedb.org/tools_api",
        "timeout": 30,
        "retries": 3,
        "retry_delay": 1
    },

    "validation": {
        "min_length": 20
    },

    "biological_data": {
        "hla_frequencies": {
            "HLA-A*02:01": 0.45,
            "HLA-A*01:01": 0.25,
            "HLA-A*03:01": 0.20,
            "HLA-B*07:02": 0.15,
            "HLA-B*08:01": 0.12,
            "HLA-C*07:01": 0.28,
            "HLA-DRB1*01:01": 0.15,
            "HLA-DRB1*15:01": 0.18
        }
    }
}