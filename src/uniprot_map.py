"""uniprot_map.py

Map ProteinGym DMS target-sequence coordinates onto UniProt canonical-sequence
coordinates so that UniProt functional annotations (active site, binding site,
domains, ...) can be transferred to DMS positions.

For each UniProt_ID (e.g. "P04637_HUMAN"):
  * accession is the part before the first underscore
  * the canonical UniProt sequence is fetched once and cached
  * target_seq (from the ProteinGym reference) is aligned to it with
    difflib.SequenceMatcher (near-identical sequences: far faster than
    Needleman-Wunsch, which is prohibitively slow for 1000+ aa proteins)
  * the resulting position map (target position -> uniprot position) is saved

Unmappable positions (target residues absent from any matching block) are
dropped and counted in the QC output.
"""
import sys
import time
import urllib.request
from pathlib import Path

import difflib

import pandas as pd

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RAW, Job  # noqa: E402


def acc_of(uniprot_id: str) -> str:
    return uniprot_id.split("_")[0]


def resolve_primary_accession(uniprot_id: str, json_cache: Path) -> str:
    """Use the primary accession from the cached UniProt JSON when available.

    Many ProteinGym UniProt_IDs are short codes (e.g. 'ACE2_HUMAN') that are
    NOT UniProt accessions; the JSON REST endpoint resolves them via search,
    and the response contains the canonical primaryAccession.
    """
    p = json_cache / f"{uniprot_id}.json"
    if p.exists():
        try:
            import json
            data = json.loads(p.read_text(encoding="utf-8"))
            pa = data.get("primaryAccession")
            if pa:
                return pa
        except Exception:
            pass
    return acc_of(uniprot_id)


def fetch_sequence(uniprot_id: str, cache: Path, job: Job, retries=5) -> str:
    p = cache / f"{uniprot_id}.fasta"
    if p.exists():
        return "".join(p.read_text().splitlines()[1:])
    acc = resolve_primary_accession(uniprot_id, cache.parent)
    url = UNIPROT_FASTA.format(acc=acc)
    for t in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "protein-ai/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                text = r.read().decode()
            lines = text.strip().splitlines()
            seq = "".join(lines[1:])
            p.write_text(text, encoding="utf-8")
            return seq
        except Exception as e:
            if t == retries - 1:
                job.info(f"uniprot fasta FAILED {uniprot_id} ({acc}): {e}")
                return ""
            time.sleep(2 * (t + 1))
    return ""


def align_map(target: str, uni_seq: str) -> dict:
    """Return {target_position (1-based): uniprot_position (1-based)}.

    Uses difflib.SequenceMatcher matching blocks; target residues outside any
    matching block are unmapped.
    """
    if not target or not uni_seq:
        return {}
    sm = difflib.SequenceMatcher(None, target, uni_seq, autojunk=False)
    out = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(i2 - i1):
                out[i1 + 1 + k] = j1 + 1 + k
    return out


def build_position_maps(target_map: dict, cache: Path, job: Job) -> pd.DataFrame:
    rows = []
    for uniprot_id, target in target_map.items():
        uni_seq = fetch_sequence(uniprot_id, cache, job)
        if not uni_seq:
            rows.append({"UniProt_ID": uniprot_id, "target_pos": None,
                         "uniprot_pos": None, "mapped": False,
                         "target_len": len(target), "uni_len": 0})
            continue
        pmap = align_map(target, uni_seq)
        for t_pos, u_pos in pmap.items():
            rows.append({"UniProt_ID": uniprot_id, "target_pos": t_pos,
                         "uniprot_pos": u_pos, "mapped": True,
                         "target_len": len(target), "uni_len": len(uni_seq)})
        job.info(f"{uniprot_id}: target {len(target)} -> uni {len(uni_seq)}, "
                 f"mapped {len(pmap)} positions "
                 f"({len(pmap) / len(target):.2f})")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    # quick self-test with a known protein
    job = Job("uniprot_map_selftest")
    cache = RAW / "uniprot" / "fasta"
    cache.mkdir(parents=True, exist_ok=True)
    seq = fetch_sequence("P04637_HUMAN", cache, job)
    print("P04637 len:", len(seq))
    m = align_map(seq[:120], seq)
    print("identity alignment first 120:", len(m), m.get(1), m.get(120))
    job.close()