# iedb_api.py

import requests
import time

def safe_post(endpoint, data, config):
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