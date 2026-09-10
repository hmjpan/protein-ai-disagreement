"""37_model_table.py

Reviewer-response (item 11-3): machine-readable model panel table.
For every core-panel model: architecture family, input modalities
(MSA / sequence / structure), whether it consumes AlphaFold structures,
and panel membership (40-model core / 3-family mechanistic).

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
ARCH_USES_STRUCT = {"ESM-IF1", "ProteinMPNN", "MIF", "MIFST", "SaProt_650M_AF2",
                    "SaProt_35M_AF2", "ProtSSN", "S2F", "S2F_MSA", "S3F",
                    "S3F_MSA", "MIFST", "ProSST", "ESCOTT", "VenusREM",
                    "RSALOR", "MULAN_small", "TranceptEVE_L", "Tranception_L",
                    "AIDO.Protein-RAG-16B"}
ARCH_USES_MSA = {"GEMME", "EVE_ensemble", "EVE_single", "DeepSequence_ensemble",
                 "DeepSequence_single", "MSA_Transformer_ensemble",
                 "MSA_Transformer_single", "EVmutation", "Site_Independent",
                 "Tranception_L", "Tranception_M", "Tranception_S",
                 "TranceptEVE_L", "TranceptEVE_M", "TranceptEVE_S",
                 "Unirep_evotune", "Wavenet", "PoET", "Protriever", "SiteRM",
                 "S2F_MSA", "S3F_MSA", "MSA-VAE", "MIFST"}


def main():
    job = Job("37_model_table")
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    models = cfg["model_list_zero_shot_substitutions_DMS"]
    panel = pd.read_csv(TABLES / "model_panel.csv")
    core = panel[panel["panel"] == "core"]["model"].tolist()
    mech = {"evolution", "single_seq", "structure"}
    rows = []
    for m in core:
        det = models.get(m, {})
        fam = det.get("model_type", "?")
        rows.append({
            "model": m,
            "official_family": fam,
            "family_in_analysis": dict(zip(panel.model, panel.family)).get(m, "?"),
            "uses_MSA": m in ARCH_USES_MSA,
            "uses_structure": (m in ARCH_USES_STRUCT
                               or "structure" in fam.lower()),
            "structure_input_is_AlphaFold": m in ARCH_USES_STRUCT,
            "in_mechanistic_panel": dict(zip(panel.model, panel.family)).get(m) in mech,
        })
    t = pd.DataFrame(rows)
    t.to_csv(STATISTICS / "model_panel_table.csv", index=False)
    job.info(f"model table written: {len(t)} models")
    job.info("structure models consuming AF input: "
             + ", ".join(t[t.structure_input_is_AlphaFold & (t.family_in_analysis == 'structure')].model))
    job.close()


if __name__ == "__main__":
    main()