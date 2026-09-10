"""Shared paths, constants, and logging for the protein AI disagreement project.

All scripts import from this module. PROJECT_ROOT can be overridden via the
environment variable to support staging/portable runs.
"""
import datetime
import os
import time
from pathlib import Path

ROOT = Path(os.environ.get(
    "PROJECT_ROOT", str(Path(__file__).resolve().parents[1])))

RAW = ROOT / "data" / "raw"
DMS_DIR = RAW / "proteingym_dms" / "DMS_ProteinGym_substitutions"
SCORES_DIR = RAW / "proteingym_scores" / "zero_shot_substitutions_scores"
CLINICAL_DIR = RAW / "proteingym_clinical" / "clinical_ProteinGym_substitutions"
STRUCTURES_DIR = RAW / "structures" / "ProteinGym_AF2_structures"
REFERENCE_DIR = RAW / "reference"

INTERMEDIATE = ROOT / "data" / "intermediate"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
STATISTICS = RESULTS / "statistics"
MODELS = RESULTS / "models"
LOGS = ROOT / "logs"
FIGURES_MAIN = ROOT / "figures" / "main"
FIGURES_SUPP = ROOT / "figures" / "supplementary"

SEED = 2026

PROTEOGYM_VERSION = "v1.3"
BASE_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3"

for d in (INTERMEDIATE, PROCESSED, RESULTS, TABLES, STATISTICS, MODELS, LOGS,
          FIGURES_MAIN, FIGURES_SUPP):
    d.mkdir(parents=True, exist_ok=True)


class Job:
    """Per-script logger writing a machine-readable run report to logs/<name>.log."""

    def __init__(self, name: str):
        self.name = name
        self.t0 = time.time()
        self.lines = []
        self.log_path = LOGS / f"{name}.log"

    def info(self, msg: str):
        line = f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        self.lines.append(line)
        print(line, flush=True)

    def section(self, title: str):
        self.info(f"=== {title} ===")

    def close(self):
        self.info(f"elapsed_seconds={time.time() - self.t0:.1f}")
        with open(self.log_path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.lines) + "\n")


def read_reference_dms() -> "pd.DataFrame":
    import pandas as pd
    return pd.read_csv(REFERENCE_DIR / "DMS_substitutions.csv")


def read_reference_clinical() -> "pd.DataFrame":
    import pandas as pd
    return pd.read_csv(REFERENCE_DIR / "clinical_substitutions.csv")


def load_core_panel(job=None) -> "pd.DataFrame":
    """Return the core model panel (panel == 'core') from model_panel.csv."""
    import pandas as pd
    p = TABLES / "model_panel.csv"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing -- run 05_quality_control.py first")
    df = pd.read_csv(p)
    core = df[df["panel"] == "core"].copy()
    if job is not None:
        job.info(f"core panel loaded: {len(core)} models")
    return core