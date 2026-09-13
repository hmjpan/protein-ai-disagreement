# -*- coding: utf-8 -*-
"""Compose main figures from panel PNGs:

1. Create the two missing Figure 1 panels: (a) workflow, (b) dataset
   composition (the only panels not produced by 16_make_figures.py).
2. Stitch Figure 1-6 panels into single labelled images
   (Fig{N}_full.png in figures/composed/), with bold lowercase letters
   matching the panel letters in the manuscript captions.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "figures" / "main"
OUTD = ROOT / "figures" / "composed"
STAT = ROOT / "results" / "statistics"
OUTD.mkdir(exist_ok=True)


def panel_file(fig, letter):
    hits = sorted(MAIN.glob(f"Fig{fig}{letter}_*.png"))
    return hits[0] if hits else None


def make_fig1a():
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    boxes = [
        (0.15, 4.6, 2.9, 1.1, "ProteinGym v1.3\n217 DMS assays\n696,311 missense variants", "#eef3fb"),
        (3.55, 4.6, 2.9, 1.1, "40-model zero-shot panel\n(evolution / sequence /\nstructure / hybrid)", "#eef3fb"),
        (6.95, 4.6, 2.9, 1.1, "Normalisation\nz-scale and percentile space\n(per assay, per model)", "#eef3fb"),
        (0.15, 2.4, 2.9, 1.4, "Model landscape\nPCA, pairwise correlations,\nregimes (GMM, K = 6)", "#eaf5ec"),
        (3.55, 2.4, 2.9, 1.4, "Structural-confidence axis\nresidue-level pLDDT gradient;\nIDR and mechanics controls", "#eaf5ec"),
        (6.95, 2.4, 2.9, 1.4, "External controls\nexperimental-structure gain;\nSSEmb replication", "#eaf5ec"),
        (0.15, 0.35, 4.6, 1.15, "Predictive probe\nBioGate CV (DMS) and\nclinical transfer (ClinVar benchmark)", "#fdf1e3"),
        (5.25, 0.35, 4.6, 1.15, "Reading\nAI-generated map of disagreement;\nexplanatory, not predictive", "#f3eefa"),
    ]
    for x, y, w, h, txt, fc in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc,
                                   edgecolor="#555555", lw=0.9))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                fontsize=8.2, linespacing=1.35)
    for x0, y0, x1, y1 in [
        (3.05, 5.15, 3.55, 5.15), (6.45, 5.15, 6.95, 5.15),
        (1.6, 4.6, 1.6, 3.8), (5.0, 4.6, 5.0, 3.8), (8.4, 4.6, 8.4, 3.8),
        (1.6, 2.4, 1.6, 1.5), (5.0, 2.4, 2.52, 1.0),
        (8.4, 2.4, 7.55, 1.0), (4.75, 0.92, 5.25, 0.92),
    ]:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.1))
    fig.tight_layout()
    fig.savefig(MAIN / "Fig1A_workflow.png", dpi=300)
    plt.close(fig)


def make_fig1b():
    mp = pd.read_csv(STAT / "model_panel_table.csv")
    fams = mp["family_in_analysis"].value_counts().reindex(
        ["evolution", "single_seq", "structure", "hybrid"]).fillna(0).astype(int)
    disp = {"evolution": "evolution", "single_seq": "single sequence",
            "structure": "structure", "hybrid": "hybrid"}
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    cols = {"evolution": "#4c72b0", "single_seq": "#dd8452",
            "structure": "#55a868", "hybrid": "#c44e52"}
    order = list(fams.index)[::-1]
    bars = ax.barh([disp[f] for f in order], list(fams.values)[::-1],
                   color=[cols[f] for f in order])
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_xlabel("core-panel models (n = 40)")
    ax.set_xlim(0, fams.max() + 3)
    ax.set_title("40-model panel by mechanistic family\n"
                 "217 assays / 186 proteins / 696,311 variants\n"
                 "(hybrid models excluded from mechanistic analyses)",
                 fontsize=10)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(MAIN / "Fig1B_dataset.png", dpi=300)
    plt.close(fig)


def get_font(size):
    for cand in ("arialbd.ttf", "segoeuib.ttf", "arial.ttf"):
        try:
            return ImageFont.truetype(cand, size)
        except OSError:
            continue
    return ImageFont.load_default()


def label(img, letter, fscale=1.0):
    im = img.convert("RGB").copy()
    d = ImageDraw.Draw(im)
    size = max(18, int(min(im.size) * 0.07 * fscale))
    font = get_font(size)
    t = f"({letter})"
    x0, y0 = 6, 4
    tb = d.textbbox((x0, y0), t, font=font)
    d.rectangle((tb[0] - 3, tb[1] - 3, tb[2] + 3, tb[3] + 3), fill="white")
    d.rectangle((tb[0] - 3, tb[1] - 3, tb[2] + 3, tb[3] + 3),
                outline="#999999", width=1)
    d.text((x0, y0), t, font=font, fill="black")
    return im


def row_height(imgs, roww):
    """Scale factor so a row of images fits width roww at shared height."""
    h = min(roww * im.size[1] / im.size[0] for im in imgs)
    return h


def compose(n, rows):
    imgs = {}
    letters = "".join("".join(r) for r in rows)
    for L in letters:
        pf = panel_file(n, L)
        if pf is None:
            raise FileNotFoundError(f"Fig{n}{L}")
        imgs[L] = Image.open(pf)
    gap = 24
    page_w = 2300
    scaled_rows = []
    for r in rows:
        ws = [imgs[L].size[0] for L in r]
        total_w = sum(ws) + gap * (len(r) - 1)
        s = (page_w - 40) / total_w
        scaled_rows.append([(imgs[L].resize((int(w * s), int(imgs[L].size[1] * s)),
                                            Image.LANCZOS), L)
                            for L, w in zip(r, ws)])
    row_h = [[im.size[1] for im, _ in sr] for sr in scaled_rows]
    W = page_w
    H = sum(max(h) for h in row_h) + gap * (len(rows) - 1) + 30
    canvas = Image.new("RGB", (W, H), "white")
    y = 15
    for sr, hs in zip(scaled_rows, row_h):
        hrow = max(hs)
        x = 20
        for im, L in sr:
            canvas.paste(im, (x, y))
            lab = label(im, L.lower())
            canvas.paste(lab, (x, y))
            x += im.size[0] + gap
        y += hrow + gap
    canvas.save(OUTD / f"Fig{n}_full.png")
    print(f"Fig{n}_full.png", canvas.size)


if __name__ == "__main__":
    for L in ("A", "B"):
        if panel_file(1, L) is None:
            (make_fig1a if L == "A" else make_fig1b)()
    compose(1, [["A", "B"], ["C", "D"], ["E"]])
    compose(2, [["A", "B"], ["C", "D"]])
    compose(3, [["A"], ["B"], ["C"]])
    compose(4, [["A"], ["B"], ["C"], ["D"]])
    compose(5, [["A"], ["B"], ["C"]])
    compose(6, [["A"], ["B"]])
    print("composed figures written to", OUTD)
