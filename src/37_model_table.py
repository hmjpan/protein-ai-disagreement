"""37_model_table.py

Machine-readable model panel table (Supplementary Table S6).
All modality metadata is DERIVED from the official ProteinGym v1.3
config.json (model_type per model), with two methodological corrections
documented in the notes column:

  - Tranception/TranceptEVE: official model_type is MSA, but the MSA is
    used only for inference-time homolog retrieval; the scored model is
    an autoregressive single-sequence transformer. Not structure-based.
  - MIF-ST: does not take MSA input; it transfers sequence representations
    from a pretrained single-sequence language model (structure is the
    generative input).

structure_input_is_AlphaFold: the released ProteinGym structure-based
scores are computed on the AlphaFold2 target structures provided by
ProteinGym (no experimental structures in the DMS benchmark).

Outputs:
  results/statistics/model_panel_table.csv
"""
import sys
import json
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, TABLES, RAW, STATISTICS, Job  # noqa: E402

CFG = RAW / "reference" / "proteingym_config.json"

# model name -> note text (methodological clarifications)
NOTES = {
    "Tranception_L": "MSA used only for inference-time homolog retrieval; "
                     "autoregressive single-sequence transformer is scored",
    "TranceptEVE_L": "MSA used only for inference-time homolog retrieval; "
                     "autoregressive transformer + EVE prior; not structure-based",
    "MIFST": "no MSA input; sequence representations transferred from a "
             "pretrained single-sequence language model",
    "ESM3": "multimodal: sequence, structure and function annotations",
    "SaProt_650M_AF2": "structure-aware vocabulary from AlphaFold2 structure",
}


def derives(model_type: str):
    uses_msa = model_type in ("MSA", "Structure & MSA")
    uses_structure = "Structure" in model_type
    return uses_msa, uses_structure


def main():
    job = Job("37_model_table")
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    models = cfg["model_list_zero_shot_substitutions_DMS"]
    panel = pd.read_csv(TABLES / "model_panel.csv")
    core = panel[panel["panel"] == "core"]["model"].tolist()
    fam_map = dict(zip(panel["model"], panel["family"]))
    mech = {"evolution", "single_seq", "structure"}
    rows = []
    for m in core:
        det = models.get(m)
        if det is None:
            job.warn(f"{m}: not in config.json")
            continue
        fam = det["model_type"]
        uses_msa, uses_structure = derives(fam)
        # Tranception: official MSA (retrieval); structure input absent.
        if m.startswith("Tranception") or m.startswith("TranceptEVE"):
            uses_structure = False
        # MIF-ST: no MSA input (transferred sequence representations only).
        if m == "MIFST":
            uses_msa = False
        rows.append({
            "model": m,
            "official_family": fam,
            "family_in_analysis": fam_map.get(m, "?"),
            "uses_MSA": uses_msa,
            "uses_structure": uses_structure,
            "structure_input_is_AlphaFold": uses_structure,
            "in_mechanistic_panel": fam_map.get(m) in mech,
            "notes": NOTES.get(m, ""),
        })
    t = pd.DataFrame(rows)
    out = STATISTICS / "model_panel_table.csv"
    t.to_csv(out, index=False)
    job.info(f"{len(t)} models -> {out.name}")
    # internal consistency audit
    for _, r in t.iterrows():
        if r["official_family"] == "Structure & MSA" and not r["uses_MSA"]:
            raise AssertionError(f"consistency: {r['model']}")
        if "Structure" not in r["official_family"] and r["uses_structure"]:
            raise AssertionError(f"structure claim without official basis: {r['model']}")



if __name__ == "__main__":
    main()
