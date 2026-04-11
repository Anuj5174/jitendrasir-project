# processing.py

import pandas as pd

from iedb_api import safe_post
from parser import parse_tsv

def get_mhci(seq, config):
    dfs = []
    for allele in config["alleles"]["mhc1"]:
        text = safe_post("mhci", {
            "method": config["prediction"]["mhc1_method"],
            "sequence_text": seq,
            "allele": allele,
            "length": config["prediction"]["mhc1_length"]
        }, config)
        df = parse_tsv(text)
        df["allele"] = allele
        dfs.append(df)

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
        text = safe_post("mhcii", {
            "method": config["prediction"]["mhc2_method"],
            "sequence_text": seq,
            "allele": allele
        }, config)
        df = parse_tsv(text)
        df["allele"] = allele
        dfs.append(df)

    df_all = pd.concat(dfs)
    if "percentile_rank" in df_all.columns:
        return df_all[df_all["percentile_rank"] <= config["thresholds"]["percentile_rank"]]
    elif "rank" in df_all.columns:
        return df_all[df_all["rank"] <= config["thresholds"]["percentile_rank"]]
    elif "ic50" in df_all.columns:
        return df_all[df_all["ic50"] < config["thresholds"]["ic50"]]
    else:
        return pd.DataFrame()


def get_bcell(seq, config):
    text = safe_post("bcell", {
        "method": config["prediction"]["bcell_method"],
        "sequence_text": seq
    }, config)
    return parse_tsv(text)