# processing.py
import pandas as pd
import requests
import time
from io import StringIO

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def validate_sequence(seq, config):
    """Clean and validate input protein sequence."""
    seq = seq.upper().strip().replace("\n", "").replace("\r", "")
    if len(seq) < config["validation"]["min_length"]:
        raise ValueError(f"Sequence too short (min {config['validation']['min_length']} AA)")
    for aa in seq:
        if aa not in VALID_AA:
            raise ValueError(f"Invalid amino acid: {aa}")
    return seq

def parse_tsv(text):
    """Simple TSV to DataFrame parser."""
    return pd.read_csv(StringIO(text), sep="\t")

def safe_post(endpoint, data, config):
    """Resilient POST request to IEDB APIs with retries."""
    base_url = config["api"]["base_url"]
    retries = config["api"]["retries"]
    timeout = config["api"]["timeout"]
    retry_delay = config["api"]["retry_delay"]
    url = f"{base_url}/{endpoint}/"

    for i in range(retries):
        try:
            res = requests.post(url, data=data, timeout=timeout)
            if res.status_code == 200 and res.text.strip():
                return res.text
        except Exception:
            time.sleep(retry_delay)
    raise Exception(f"IEDB API failed: {endpoint}")

def get_mhci(seq, config):
    # (Rest of get_mhci remains same)
    # ...
    dfs = []
    for allele in config["alleles"]["mhc1"]:
        try:
            text = safe_post("mhci", {
                "method": config["prediction"]["mhc1_method"],
                "sequence_text": seq,
                "allele": allele,
                "length": config["prediction"]["mhc1_length"]
            }, config)
            df = parse_tsv(text)
            df["allele"] = allele
            dfs.append(df)
        except Exception as e:
            print(f"Warning: MHC-I prediction failed for allele {allele}: {e}")
            continue

    df_all = pd.concat(dfs)
    if "percentile_rank" in df_all.columns:
        return df_all[df_all["percentile_rank"] <= config["thresholds"]["percentile_rank"]]
    elif "rank" in df_all.columns:
        return df_all[df_all["rank"] <= config["thresholds"]["percentile_rank"]]
    elif "ic50" in df_all.columns:
        return df_all[df_all["ic50"] < config["thresholds"]["ic50"]]
    else:
        return pd.DataFrame()


def get_mhcii(seq, config):
    dfs = []
    for allele in config["alleles"]["mhc2"]:
        try:
            text = safe_post("mhcii", {
                "method": config["prediction"]["mhc2_method"],
                "sequence_text": seq,
                "allele": allele
            }, config)
            df = parse_tsv(text)
            df["allele"] = allele
            dfs.append(df)
        except Exception as e:
            print(f"Warning: MHC-II prediction failed for allele {allele}: {e}")
            continue

    df_all = pd.concat(dfs)
    if "percentile_rank" in df_all.columns:
        return df_all[df_all["percentile_rank"] <= config["thresholds"]["percentile_rank"]]
    elif "rank" in df_all.columns:
        return df_all[df_all["rank"] <= config["thresholds"]["percentile_rank"]]
    elif "ic50" in df_all.columns:
        ic50_thresh = config["thresholds"].get("mhc2_ic50", config["thresholds"]["ic50"])
        return df_all[df_all["ic50"] < ic50_thresh]
    else:
        return pd.DataFrame()


def get_bcell(seq, config):
    try:
        text = safe_post("bcell", {
            "method": config["prediction"]["bcell_method"],
            "sequence_text": seq
        }, config)
        return parse_tsv(text)
    except Exception as e:
        print(f"Warning: B-cell prediction failed: {e}")
        return pd.DataFrame()