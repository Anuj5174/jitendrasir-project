import sys
import os

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

import uvicorn
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles


# Add module directories to sys.path so imports work correctly
MODULE1_DIR = os.path.join(os.path.dirname(__file__), "module1")
MODULE2_DIR = os.path.join(os.path.dirname(__file__), "module2")

sys.path.insert(0, MODULE1_DIR)
sys.path.insert(0, MODULE2_DIR)

# Import the FastAPI apps from both modules
from module1.module1_api import app as module1_app
from module2.structure_api import app as module2_app

app = FastAPI()

# Mount module2 on /module2 to expose its API (e.g. /module2/api/fold)
# and its static assets.
app.mount("/module2", module2_app)

# Mount module1 on the root level so the main dashboard is accessible at /
app.mount("/", module1_app)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8080))
    uvicorn.run(app, host="0.0.0.0", port=port)
