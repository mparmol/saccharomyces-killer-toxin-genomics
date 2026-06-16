#!/usr/bin/env python3
"""
16c_colored_leaves.py
NJ unrooted tree — hojas coloreadas por especie, sin cajas.
nhmmer candidates marcados con * en rojo/itálica.
Output: results_2/networks_colored/
"""

import re
import numpy as np
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import dendropy

PROJECT   = Path("")
RESULTS   = PROJECT / "results_2"
OUT_DIR   = RESULTS / "networks_colored"
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

def nj_layout(taxa, D):
    n = len(taxa)
    tns = dendropy.TaxonNamespace(taxa)
    tl  = [tns.get_taxon(s) for s in taxa]
    dd  = {tl[i]: {tl[j]: D[i, j] for j in range(n) if i != j} for i in range(n)}
    pdm = dendropy.PhylogeneticDistanceMatrix()
    pdm.compile_from_dict(dd, taxon_namespace=tns)
    tree = pdm.nj_tree()
    tree.reroot_at_midpoint(update_bipartitions=True)
    raw = [nd.edge_length for nd in tree.preorder_node_iter()
           if nd.edge_length is not None and nd.edge_length > 0]
    if raw:
        med = np.median(raw)
        cap = med * 6
        def compress(l):
            l = max(abs(l or 1e-6), 1e-6)
            return np.sqrt(min(l, cap) / (med + 1e-9)) * med
    else:
        compress = lambda l: 0.01
    for nd in tree.preorder_node_iter():
        nd.edge_length = max(compress(nd.edge_length or 1e-6), 1e-4)
    lc = {}
    def count_leaves(nd):
        if nd.is_leaf():
            lc[nd] = 1
        else:
            for c in nd.child_nodes():
                count_leaves(c)
            lc[nd] = sum(lc[c] for c in nd.child_nodes())
    count_leaves(tree.seed_node)
    pos = {}
    def place(nd, start_angle, arc, parent_xy):
        mid = start_angle + arc / 2
        el  = nd.edge_length
        x   = parent_xy[0] + el * np.cos(mid)
        y   = parent_xy[1] + el * np.sin(mid)
        pos[nd] = np.array([x, y])
        if nd.is_leaf():
            return
        children = list(nd.child_nodes())
        total    = lc[nd]
        cur      = start_angle
        for ch in children:
            ch_arc = arc * lc[ch] / total
            place(ch, cur, ch_arc, (x, y))
            cur += ch_arc
    root = tree.seed_node
    pos[root] = np.zeros(2)
    children = list(root.child_nodes())
    total = lc[root]
    cur = 0.0
    for ch in children:
        place(ch, cur, 2 * np.pi * lc[ch] / total, (0.0, 0.0))
        cur += 2 * np.pi * lc[ch] / total
    tip_rs = [np.linalg.norm(pos[nd]) for nd in tree.leaf_node_iter() if nd in pos]
    scale  = np.percentile(tip_rs, 90) if tip_rs else 1.0
    if scale < 1e-9:
        scale = 1.0
    pos = {nd: p / scale for nd, p in pos.items()}
    edges = []
    for nd in tree.preorder_node_iter():
        if nd.parent_node is not None and nd in pos and nd.parent_node in pos:
            p0, p1 = pos[nd.parent_node], pos[nd]
            edges.append((p0[0], p0[1], p1[0], p1[1]))
    leaf_pos = {nd.taxon.label: pos[nd] for nd in tree.leaf_node_iter() if nd in pos}
    return leaf_pos, edges


def draw_network(gene_label, mldist_path, lookup, species2hex,
                 detect_methods=None, show_detect=False):
    if not mldist_path.exists():
        print(f"  [SKIP] {mldist_path.name} not found")
        return

    taxa, D = parse_mldist(mldist_path)
    print(f"  [{gene_label}] {len(taxa)} taxa")
    pos, edges = nj_layout(taxa, D)

    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color="#aaaaaa", linewidth=0.55, zorder=1)

    species_seen = {}
    for t in taxa:
        if t not in pos:
            continue
        p  = pos[t]
        sp = lookup(base_name(t))
        col = species2hex.get(sp, "#888888")

        method   = detect_methods.get(base_name(t), "BLASTn") if detect_methods else "BLASTn"
        is_novel = show_detect and method == "nhmmer"

        marker_size = 28 if is_novel else 18
        edge_col    = "#cc0000" if is_novel else "white"
        ax.scatter(p[0], p[1], c=col, s=marker_size,
                   edgecolors=edge_col, linewidths=0.8 if is_novel else 0.3,
                   zorder=4, marker="o")

        all_x = np.array([pos[tt][0] for tt in taxa if tt in pos])
        all_y = np.array([pos[tt][1] for tt in taxa if tt in pos])
        cx, cy = all_x.mean(), all_y.mean()
        dx, dy = p[0] - cx, p[1] - cy
        dist   = np.sqrt(dx**2 + dy**2)
        if dist > 1e-9:
            nx, ny = dx / dist, dy / dist
        else:
            nx, ny = 1.0, 0.0
        offset = 0.04
        lx, ly = p[0] + nx * offset, p[1] + ny * offset
        ha = "left" if nx > 0 else ("right" if nx < -0.1 else "center")
        va = "bottom" if ny > 0.3 else ("top" if ny < -0.3 else "center")

        label = f"{base_name(t)}*" if is_novel else base_name(t)
        ax.text(lx, ly, label,
                ha=ha, va=va, fontsize=5.0,
                fontstyle="italic" if is_novel else "normal",
                color="#cc0000" if is_novel else "#222222",
                zorder=5)

        if sp not in species_seen:
            species_seen[sp] = col

    spp_patches = [
        mpatches.Patch(color=col, label=f"S. {sp.split()[1]}" if len(sp.split()) >= 2 else sp)
        for sp, col in sorted(species_seen.items())
        if sp != "Unknown"
    ]
    ax.legend(handles=spp_patches, loc="lower left", fontsize=8,
              framealpha=0.92, edgecolor="#cccccc", ncol=2,
              title="Especie", title_fontsize=8)

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

    ax.text(0.02, 0.98, gene_label.replace("_extended", ""),
            transform=ax.transAxes, fontsize=20, fontweight="bold",
            va="top", ha="left")
    if "extended" in gene_label.lower():
        ax.text(0.02, 0.93, "(+ nhmmer candidates)",
                transform=ax.transAxes, fontsize=9, va="top", ha="left", color="#666")

    all_x = np.array([pos[t][0] for t in taxa if t in pos])
    all_y = np.array([pos[t][1] for t in taxa if t in pos])
    cx_all, cy_all = all_x.mean(), all_y.mean()
    tip_r = np.sqrt((all_x - cx_all)**2 + (all_y - cy_all)**2)
    margin = np.median(tip_r) * 0.60
    ax.set_xlim(all_x.min() - margin, all_x.max() + margin)
    ax.set_ylim(all_y.min() - margin, all_y.max() + margin)

    plt.tight_layout(pad=0.3)
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"{gene_label}_colored.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else 150,
                    bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.close(fig)


def main():
    print("[16c] Loading metadata...")
    lookup, species2hex = load_metadata()
    detect_methods = load_detect_methods()
    for gene_label, mldist_path in DATASETS.items():
        print(f"\n[16c] {gene_label}")
        draw_network(gene_label, mldist_path, lookup, species2hex,
                     detect_methods=detect_methods,
                     show_detect="extended" in gene_label.lower())
    print("\n[16c] Done. Output:", OUT_DIR)

if __name__ == "__main__":
    main()
