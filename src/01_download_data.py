"""01_download_data.py

Verify integrity of the downloaded ProteinGym v1.3 data and record a manifest
(checksums). Fetch the official `config.json` from the ProteinGym GitHub
repository, which is the authoritative source for:
  * which models are in the zero-shot DMS substitution benchmark
  * each model's `directionality` (how raw scores were flipped so that the
    RELEASED merged score files have: higher score = higher fitness)
  * each model's `model_type` (MSA / Single sequence / Structure / hybrids)
    used as the official model-family label.

Direction is therefore taken from the official configuration only -- never
inferred from DMS labels (pre-registered rule, protocol section 10/46).

Inputs : data/raw/** zips + reference csvs (must already exist)
         https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/config.json
Outputs: data/raw/MANIFEST.md
         data/processed/model_metadata.csv   (official per-model info)
         data/processed/model_direction.csv  (direction for normalization)
"""
import hashlib
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    CLINICAL_DIR, DMS_DIR, PROCESSED, RAW, REFERENCE_DIR, SCORES_DIR,
    STRUCTURES_DIR, Job,
)

MANIFEST = [
    ("DMS_ProteinGym_substitutions.zip", "proteingym_dms", 218,
     "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip"),
    ("zero_shot_substitutions_scores.zip", "proteingym_scores", 217,
     "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip"),
    ("clinical_ProteinGym_substitutions.zip", "proteingym_clinical", 2525,
     "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/clinical_ProteinGym_substitutions.zip"),
    ("ProteinGym_AF2_structures.zip", "structures", 199,
     "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/ProteinGym_AF2_structures.zip"),
]

CONFIG_URL = ("https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/"
              "config.json")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_config(job: Job) -> dict:
    try:
        req = urllib.request.Request(CONFIG_URL, headers={"User-Agent": "protein-ai"})
        with urllib.request.urlopen(req, timeout=90) as r:
            cfg = json.loads(r.read().decode())
        job.info(f"fetched official config.json ({len(cfg)} top-level keys)")
        return cfg
    except Exception as e:
        job.info(f"config.json fetch FAILED: {e}")
        return {}


def main():
    job = Job("01_download_data")
    job.section("manifest checks")
    rows = []
    ok = True
    for fname, subdir, expected_entries, url in MANIFEST:
        p = RAW / subdir / fname
        if not p.exists():
            job.info(f"MISSING {p}")
            ok = False
            rows.append({"file": fname, "path": str(p), "url": url, "entries": -1,
                         "expected_entries": expected_entries, "status": "MISSING",
                         "sha256": ""})
            continue
        try:
            with zipfile.ZipFile(p) as z:
                n = len(z.namelist())
                bad = z.testzip()
            status = "OK" if (n == expected_entries and bad is None) else f"CORRUPT({bad})"
            if bad is not None:
                ok = False
        except zipfile.BadZipFile as e:
            n, status = -1, f"CORRUPT({e})"
            ok = False
        job.info(f"{fname}: entries={n} expected={expected_entries} status={status}")
        rows.append({"file": fname, "path": str(p), "url": url, "entries": n,
                     "expected_entries": expected_entries, "status": status,
                     "sha256": sha256(p)})

    for fname in ("DMS_substitutions.csv", "clinical_substitutions.csv"):
        p = REFERENCE_DIR / fname
        exists = p.exists()
        job.info(f"reference {fname}: exists={exists}")
        rows.append({"file": fname, "path": str(p), "url": "", "entries": -1,
                     "expected_entries": -1, "status": "OK" if exists else "MISSING",
                     "sha256": sha256(p) if exists else ""})
        if not exists:
            ok = False

    for d, label in [(DMS_DIR, "DMS"), (SCORES_DIR, "scores"),
                     (CLINICAL_DIR, "clinical"), (STRUCTURES_DIR, "structures")]:
        n = len(list(d.glob("*"))) if d.exists() else -1
        job.info(f"extracted dir {label}: {n} items")
        if n < 1:
            ok = False

    pd.DataFrame(rows).to_csv(RAW / "MANIFEST.md", index=False)

    job.section("official model metadata")
    cfg = fetch_config(job)
    models = cfg.get("model_list_zero_shot_substitutions_DMS", {})
    meta_rows = []
    dir_rows = []
    for name, det in models.items():
        directionality = det.get("directionality")
        model_type = det.get("model_type")
        # Released (merged) score files already carry raw * directionality,
        # i.e. higher score = higher fitness for every model (official merge.py).
        family = {"MSA": "evolution", "Single sequence": "single_seq",
                  "Structure": "structure"}.get(model_type, "hybrid")
        meta_rows.append({"model_name": name, "model_type": model_type,
                          "family": family, "directionality": directionality,
                          "input_score_name": det.get("input_score_name")})
        dir_rows.append({"model_name": name,
                         "released_higher_means": "fitness",
                         "mutational_direction": "higher_is_better",
                         "raw_directionality": directionality,
                         "evidence": "official_config.json"})
    if not meta_rows:
        job.info("WARNING: no model metadata fetched; pipeline cannot continue")
        ok = False
    pd.DataFrame(meta_rows).to_csv(PROCESSED / "model_metadata.csv", index=False)
    pd.DataFrame(dir_rows).to_csv(PROCESSED / "model_direction.csv", index=False)
    job.info(f"model_metadata.csv: {len(meta_rows)} models; "
             f"families: {pd.DataFrame(meta_rows).family.value_counts().to_dict()}")
    job.info(f"ALL_OK={ok}")
    job.close()


if __name__ == "__main__":
    main()