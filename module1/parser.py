# parser.py

import pandas as pd
from io import StringIO

def parse_tsv(text):
    return pd.read_csv(StringIO(text), sep="\t")