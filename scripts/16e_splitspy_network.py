#!/usr/bin/env python3
"""
16e_splitspy_network.py
Proper NeighborNet phylogenetic outline using splitspy (Bryant & Moulton 2004).
Species bounding boxes + colored leaf labels + nhmmer asterisks.
Output: results_2/networks_outline/
"""

import re
import math
import numpy as np
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from splitspy.nnet.nnet_algo import neighbor_net
import splitspy.outlines.outline_algo as outline_algo

PROJECT   = Path("")
RESULTS   = PROJECT / "results_2"
OUT_DIR   = RESULTS / "networks_outline"
OUT_DIR.mkdir(parents=True, exist_ok=True)

META_TSV   = PROJECT / "results" / "genome_metadata.tsv"
SPP_CSV    = PROJECT / "Spp_Color_codes.csv"
DETECT_TSV = RESULTS / "discovery" / "detection_method_summary.tsv"

DATASETS = {
    "KHS1":          PROJECT / "results" / "trees" / "KHS1" / "KHS1_primary_iqtree.mldist",
    "KHR1":          PROJECT / "results" / "trees" / "KHR1" / "KHR1_iqtree.mldist",
    "KHR1_extended": RESULTS / "discovery" / "combined" / "KHR1_extended_iqtree.mldist",
}

METHOD_PRIORITY = {"nhmmer": 3, "tBLASTx": 2, "BLASTn": 1}


def norm_strain(s):
    return re.sub(r'[-_.\s]', '', s).upper()

def base_name(name):
    return name.split("__")[0] if "__" in name else name

def load_metadata():
    meta = pd.read_csv(META_TSV, sep="\t")
    spp  = pd.read_csv(SPP_CSV)
    strain2species, norm2species = {}, {}
    for _, row in meta.iterrows():
        s  = str(row["Strain"]).strip()
        sp = str(row["Species"]).strip()
        strain2species[s] = sp
        norm2species[norm_strain(s)] = sp
    species2hex = {}
    for _, row in spp.iterrows():
        sp = str(row["Species"]).strip()
        species2hex[sp] = str(row["HEX code"]).strip()
    def lookup(name):
        name = base_name(name)
        if name in strain2species:
            return strain2species[name]
        n = norm_strain(name)
        if n in norm2species:
            return norm2species[n]
        if n and n[-1].isalpha() and n[:-1] in norm2species:
            return norm2species[n[:-1]]
        return "Unknown"
    return lookup, species2hex

def load_detect_methods():
    if not DETECT_TSV.exists():
        return {}
    dt = pd.read_csv(DETECT_TSV, sep="\t")
    det = {}
    for _, row in dt.iterrows():
        s, m = str(row["strain"]).strip(), str(row["detection_method"]).strip()
        if METHOD_PRIORITY.get(m, 0) > METHOD_PRIORITY.get(det.get(s, ""), 0):
            det[s] = m
    return det

def parse_mldist(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    taxa, rows = [], []
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        taxa.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    D = np.array(rows, dtype=float)
    np.fill_diagonal(D, 0.0)
    return taxa, D

def hex_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def build_splitspy_graph(taxa, D):
    # Sqrt compression: reduces dominance of long-branch outliers
    D_c = np.sqrt(D)
    np.fill_diagonal(D_c, 0.0)
    mat = D_c.tolist()

    cycle, splits = neighbor_net(taxa, mat)
    graph, angles = outline_algo.compute(taxa, cycle, splits)
    return graph, angles

def extract_graph_data(graph):
    """Returns edges as list of (x0,y0,x1,y1) and leaf dict {label: [x,y]}."""
    edges = []
    for e in graph.edges():
        p0 = e.src().pos
        p1 = e.tar().pos
        edges.append((float(p0[0]), float(p0[1]), float(p1[0]), float(p1[1])))

    leaf_pos = {}
    for v in graph.nodes():
        if v.label is not None and v.label != "Root":
            p = v.pos
            # label can be comma-separated if multiple taxa overlap (rare)
            for taxon in v.label.split(","):
                taxon = taxon.strip()
                if taxon:
                    leaf_pos[taxon] = [float(p[0]), float(p[1])]
    return edges, leaf_pos


def draw_network(gene_label, mldist_path, lookup, species2hex,
                 detect_methods=None, show_detect=False):
    if not mldist_path.exists():
        print(f"  [SKIP] {mldist_path.name} not found")
        return

    taxa, D = parse_mldist(mldist_path)
    print(f"  [{gene_label}] {len(taxa)} taxa — running NeighborNet outline...")

    graph, angles = build_splitspy_graph(taxa, D)
    edges, leaf_pos = extract_graph_data(graph)
    print(f"  [{gene_label}] {len(edges)} edges, {len(leaf_pos)} leaves in outline")

    fig, ax = plt.subplots(figsize=(16, 16))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── backbone edges ─────────────────────────────────────────────────────────
    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color="#aaaaaa", linewidth=0.55, zorder=1)

    # ── leaves: colored dots + radial labels ───────────────────────────────────
    all_x = np.array([p[0] for p in leaf_pos.values()])
    all_y = np.array([p[1] for p in leaf_pos.values()])
    cx, cy = all_x.mean(), all_y.mean()

    species_seen = {}
    for t in taxa:
        bn = base_name(t)
        if bn not in leaf_pos:
            continue
        p = leaf_pos[bn]
        sp  = lookup(bn)
        col = species2hex.get(sp, "#888888")

        method   = detect_methods.get(bn, "BLASTn") if detect_methods else "BLASTn"
        is_novel = show_detect and method == "nhmmer"

        ms = 32 if is_novel else 20
        ec = "#cc0000" if is_novel else "white"
        ax.scatter(p[0], p[1], c=col, s=ms,
                   edgecolors=ec, linewidths=0.9 if is_novel else 0.3,
                   zorder=4, marker="o")

        # radial label direction
        dx, dy = p[0] - cx, p[1] - cy
        dist = math.hypot(dx, dy)
        nx, ny = (dx / dist, dy / dist) if dist > 1e-9 else (1.0, 0.0)
        offset = (all_x.max() - all_x.min()) * 0.025
        lx, ly = p[0] + nx * offset, p[1] + ny * offset
        ha = "left" if nx > 0.1 else ("right" if nx < -0.1 else "center")
        va = "bottom" if ny > 0.3 else ("top" if ny < -0.3 else "center")

        label = f"{bn}*" if is_novel else bn
        ax.text(lx, ly, label, ha=ha, va=va, fontsize=4.5,
                fontstyle="italic" if is_novel else "normal",
                color="#cc0000" if is_novel else "#222222",
                zorder=5)

        if sp not in species_seen:
            species_seen[sp] = col

    # ── species legend (bottom-left) ───────────────────────────────────────────
    spp_patches = [
        mpatches.Patch(color=col, label=f"S. {sp.split()[1]}" if len(sp.split()) >= 2 else sp)
        for sp, col in sorted(species_seen.items())
        if sp != "Unknown"
    ]
    leg1 = ax.legend(handles=spp_patches, loc="lower left", fontsize=8,
                     framealpha=0.92, edgecolor="#cccccc", ncol=2,
                     title="Especie", title_fontsize=8)
    ax.add_artist(leg1)

    # ── detection legend (bottom-right) ───────────────────────────────────────
    if show_detect and detect_methods:
        methods_seen = set(detect_methods.get(base_name(t), "BLASTn") for t in taxa)
        if "nhmmer" in methods_seen:
            det_patches = [
                mpatches.Patch(facecolor="white", edgecolor="#aaaaaa", label="BLAST (original)"),
                mpatches.Patch(facecolor="white", edgecolor="#cc0000", label="nhmmer (novel)*"),
            ]
            ax.legend(handles=det_patches, loc="lower right", fontsize=8,
                      framealpha=0.92, edgecolor="#cccccc",
                      title="Detección", title_fontsize=8)

    # ── title ─────────────────────────────────────────────────────────────────
    ax.text(0.02, 0.98, gene_label.replace("_extended", ""),
            transform=ax.transAxes, fontsize=20, fontweight="bold",
            va="top", ha="left")
    if "extended" in gene_label.lower():
        ax.text(0.02, 0.93, "(+ nhmmer candidates)",
                transform=ax.transAxes, fontsize=9, va="top", ha="left", color="#666")
    ax.text(0.02, 0.02, "NeighborNet outline (Bryant & Moulton 2004)",
            transform=ax.transAxes, fontsize=7, va="bottom", ha="left", color="#888")

    # ── axis limits ───────────────────────────────────────────────────────────
    all_xe = np.array([x for x0,y0,x1,y1 in edges for x in (x0, x1)])
    all_ye = np.array([y for x0,y0,x1,y1 in edges for y in (y0, y1)])
    xpad = (all_xe.max() - all_xe.min()) * 0.18
    ypad = (all_ye.max() - all_ye.min()) * 0.18
    ax.set_xlim(all_xe.min() - xpad, all_xe.max() + xpad)
    ax.set_ylim(all_ye.min() - ypad, all_ye.max() + ypad)

    plt.tight_layout(pad=0.3)
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"{gene_label}_outline.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else 150,
                    bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.close(fig)


def main():
    print("[16e] Loading metadata...")
    lookup, species2hex = load_metadata()
    detect_methods = load_detect_methods()
    for gene_label, mldist_path in DATASETS.items():
        print(f"\n[16e] {gene_label}")
        draw_network(gene_label, mldist_path, lookup, species2hex,
                     detect_methods=detect_methods,
                     show_detect="extended" in gene_label.lower())
    print("\n[16e] Done. Output:", OUT_DIR)

if __name__ == "__main__":
    main()
