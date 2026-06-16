#!/usr/bin/env python3
"""
17_networks_v2.py
Three network styles per dataset, with optional style config.

Styles:
  outline   — NeighborNet (splitspy), coloured tips
  nj_boxes  — NJ unrooted tree, species bounding boxes
  colored   — NJ unrooted tree, coloured tips + legend

All functions accept cfg=None (hardcoded defaults) or a dict loaded from
style_config.json for publication-quality output.  SVG is always saved
(individually selectable elements in Inkscape).

Output: results_2/networks_{outline,nj,colored}/
"""

import re
import json
import math
import numpy as np
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import dendropy

from splitspy.nnet.nnet_algo import neighbor_net
import splitspy.outlines.outline_algo as outline_algo

PROJECT   = Path("")
RESULTS   = PROJECT / ""

META_TSV   = RESULTS / "01_metadata" / "genome_metadata.tsv"
SPP_CSV    = PROJECT / "Spp_Color_codes.csv"
DETECT_TSV = RESULTS / "05_hmm_discovery" / "detection_method_summary.tsv"
STYLE_CONFIG_PATH = Path(__file__).parent / "style_config.json"

DATASETS = {
    "KHS1":          (RESULTS / "" / "KHS1" / "KHS1_primary_iqtree.mldist", "KHS1"),
    "KHR1":          (RESULTS / "" / "KHR1" / "KHR1_iqtree.mldist",         "KHR1"),
    "KHR1_extended": (RESULTS / "05_hmm_discovery" / "candidatos_KHR1" / "KHR1_extended_iqtree.mldist", "KHR1"),
}

METHOD_PRIORITY = {"nhmmer": 3, "tBLASTx": 2, "BLASTn": 1}

STRAIN_CORRECTIONS = {
    "UFRJ50816T": "UFRJ50816",
    "EM14S01-3B": "EM14S013B",
}


# ─── Data loaders ─────────────────────────────────────────────────────────────

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
        if s in STRAIN_CORRECTIONS:
            corrected = STRAIN_CORRECTIONS[s]
            strain2species[corrected] = sp
            norm2species[norm_strain(corrected)] = sp
    species2hex = {}
    for _, row in spp.iterrows():
        sp = str(row["Species"]).strip()
        species2hex[sp] = str(row["HEX code"]).strip()
    def lookup(name):
        bn = base_name(name)
        if bn in strain2species:
            return strain2species[bn]
        n = norm_strain(bn)
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

def load_chrom_lookup():
    if not DETECT_TSV.exists():
        return {}
    dt = pd.read_csv(DETECT_TSV, sep="\t")
    cl = {}
    for _, row in dt.iterrows():
        s = str(row["strain"]).strip()
        g = str(row["gene"]).strip()
        c = str(row["chrom"]).strip()
        k = (s, g)
        if k not in cl:
            cl[k] = c
    return cl

def load_style_config():
    with open(STYLE_CONFIG_PATH) as f:
        return json.load(f)

def parse_mldist(path):
    lines = Path(path).read_text().splitlines()
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


# ─── Style param resolution ────────────────────────────────────────────────────

def _p(cfg, *keys, default):
    """Safely traverse cfg dict; return default if cfg is None or key missing."""
    if cfg is None:
        return default
    val = cfg
    for k in keys:
        if not isinstance(val, dict) or k not in val:
            return default
        val = val[k]
    return val


# ─── Label builder ─────────────────────────────────────────────────────────────

def build_label(t, gene, chrom_lookup, detect_methods):
    if '__' in t:
        strain    = t.split('__')[0]
        chrom     = t.split('__')[1].split('_')[0]
        is_nhmmer = 'nhmmer' in t
    else:
        strain    = t
        chrom     = chrom_lookup.get((strain, gene), "")
        is_nhmmer = False
    chrom_s = chrom.replace('chr', '') if chrom else '?'
    mark    = '*' if is_nhmmer else ''
    return f"{strain}{mark} ({chrom_s})", strain, is_nhmmer


def radial_label_pos(px, py, cx, cy, xrange):
    dx, dy = px - cx, py - cy
    dist = math.hypot(dx, dy)
    nx, ny = (dx / dist, dy / dist) if dist > 1e-9 else (1.0, 0.0)
    offset = xrange * 0.025
    lx, ly = px + nx * offset, py + ny * offset
    angle = math.degrees(math.atan2(dy, dx))
    if nx < -0.1:      # left half — flip so text reads left-to-right
        angle += 180
        ha = "right"
    elif nx > 0.1:
        ha = "left"
    else:
        ha = "center"
    return lx, ly, ha, "center", angle


# ─── NJ layout ────────────────────────────────────────────────────────────────

def nj_layout(taxa, D):
    n   = len(taxa)
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


# ─── splitspy outline ─────────────────────────────────────────────────────────

def build_splitspy_graph(taxa, D):
    D_c = np.sqrt(D)
    np.fill_diagonal(D_c, 0.0)
    cycle, splits = neighbor_net(taxa, D_c.tolist())
    graph, angles = outline_algo.compute(taxa, cycle, splits)
    return graph, angles

def extract_graph_data(graph):
    edges = []
    for e in graph.edges():
        p0 = e.src().pos
        p1 = e.tar().pos
        edges.append((float(p0[0]), float(p0[1]), float(p1[0]), float(p1[1])))
    leaf_pos = {}
    for v in graph.nodes():
        if v.label is not None and v.label != "Root":
            p = v.pos
            for taxon in v.label.split(","):
                taxon = taxon.strip()
                if taxon:
                    leaf_pos[taxon] = [float(p[0]), float(p[1])]
    return edges, leaf_pos


# ─── Shared legend / title helpers ────────────────────────────────────────────

def _species_legend(ax, species_seen, loc="lower left", cfg=None):
    fs  = _p(cfg, "legend", "fontsize",       default=9)
    tfs = _p(cfg, "legend", "title_fontsize", default=9)
    nc  = _p(cfg, "legend", "ncol",           default=2)
    patches = [
        mpatches.Patch(color=col, label=f"S. {sp.split()[1]}" if len(sp.split()) >= 2 else sp)
        for sp, col in sorted(species_seen.items()) if sp != "Unknown"
    ]
    return ax.legend(handles=patches, loc=loc, framealpha=0.92, edgecolor="#cccccc",
                     ncol=nc, title="Species", title_fontsize=tfs,
                     prop={"size": fs, "style": "italic"})

def _detection_legend(ax, loc="lower right", cfg=None):
    fs  = _p(cfg, "legend", "fontsize",       default=9)
    tfs = _p(cfg, "legend", "title_fontsize", default=9)
    patches = [
        mpatches.Patch(facecolor="white", edgecolor="#aaaaaa", label="BLAST (original)"),
        mpatches.Patch(facecolor="white", edgecolor="#cc0000", label="nhmmer (novel)*"),
    ]
    ax.legend(handles=patches, loc=loc, fontsize=fs,
              framealpha=0.92, edgecolor="#cccccc",
              title="Detection", title_fontsize=tfs)

def _draw_title(ax, gene_label, cfg=None):
    fs = _p(cfg, "title", "fontsize",   default=22)
    fw = _p(cfg, "title", "fontweight", default="bold")
    ax.text(0.02, 0.98, gene_label.replace("_extended", ""),
            transform=ax.transAxes, fontsize=fs, fontweight=fw,
            fontstyle="italic", va="top", ha="left")
    if "extended" in gene_label.lower():
        ax.text(0.02, 0.94, "(+ nhmmer candidates)",
                transform=ax.transAxes, fontsize=max(fs * 0.25, 10),
                va="top", ha="left", color="#666")

def _has_nhmmer(taxa):
    return any("nhmmer" in t for t in taxa)

def _save(fig, out_dir, stem):
    for ext in ("pdf", "png", "svg"):
        out = out_dir / f"{stem}.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else 150,
                    bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")


# ─── Style 1: NeighborNet outline ─────────────────────────────────────────────

def draw_outline(gene_label, gene, mldist_path, out_dir,
                 lookup, species2hex, chrom_lookup, detect_methods,
                 cfg=None, suffix=""):
    if not mldist_path.exists():
        print(f"  [SKIP] {mldist_path.name}")
        return

    taxa, D = parse_mldist(mldist_path)
    print(f"  [outline/{gene_label}] {len(taxa)} taxa — running NeighborNet...")
    graph, _ = build_splitspy_graph(taxa, D)
    edges, leaf_pos = extract_graph_data(graph)
    print(f"  [outline/{gene_label}] {len(edges)} edges, {len(leaf_pos)} leaves")

    figsize    = tuple(_p(cfg, "figure", "figsize", default=(20, 20)))
    edge_col   = _p(cfg, "edges", "color",     default="#aaaaaa")
    edge_lw    = _p(cfg, "edges", "linewidth", default=0.6)
    s_blast    = _p(cfg, "markers", "size_blast",       default=55)
    s_nhmmer   = _p(cfg, "markers", "size_nhmmer",      default=80)
    ec_nhmmer  = _p(cfg, "markers", "edge_color_nhmmer",default="#cc0000")
    lw_blast   = _p(cfg, "markers", "linewidth_blast",  default=0.4)
    lw_nhmmer  = _p(cfg, "markers", "linewidth_nhmmer", default=1.0)
    lbl_fs     = _p(cfg, "labels", "fontsize",    default=7.0)
    lbl_col_nm = _p(cfg, "labels", "color_nhmmer",default="#cc0000")
    highlights = _p(cfg, "highlights", default={}) if cfg else {}
    show_meth  = _p(cfg, "method_text", "show", default=True)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color=edge_col, linewidth=edge_lw, zorder=1)

    all_lx = np.array([p[0] for p in leaf_pos.values()])
    all_ly = np.array([p[1] for p in leaf_pos.values()])
    cx, cy = all_lx.mean(), all_ly.mean()
    xrange = all_lx.max() - all_lx.min() if len(all_lx) else 1.0

    species_seen, missing = {}, []
    for t in taxa:
        if t not in leaf_pos:
            missing.append(t)
            continue
        p  = leaf_pos[t]
        sp = lookup(t)
        col = species2hex.get(sp, "#888888")
        label, strain, is_nhmmer = build_label(t, gene, chrom_lookup, detect_methods)

        hl       = highlights.get(strain, {})
        marker_s = hl.get("marker_size", s_nhmmer if is_nhmmer else s_blast)
        lbl_size = hl.get("label_fontsize", lbl_fs)
        lbl_fw   = hl.get("label_fontweight", "normal")

        ax.scatter(p[0], p[1], c=col, s=marker_s,
                   edgecolors=ec_nhmmer if is_nhmmer else "white",
                   linewidths=lw_nhmmer if is_nhmmer else lw_blast,
                   zorder=4, marker="o")

        lx, ly, ha, va, angle = radial_label_pos(p[0], p[1], cx, cy, xrange)
        ax.text(lx, ly, label, ha=ha, va=va,
                fontsize=lbl_size, fontweight=lbl_fw,
                fontstyle="italic" if is_nhmmer else "normal",
                color=lbl_col_nm if is_nhmmer else col,
                rotation=angle, rotation_mode="anchor",
                zorder=5)

        if sp not in species_seen:
            species_seen[sp] = col

    if missing:
        print(f"  WARNING — {len(missing)} missing from leaf_pos: {missing[:8]}")

    leg1 = _species_legend(ax, species_seen, "lower left", cfg)
    ax.add_artist(leg1)
    if "extended" in gene_label.lower() and _has_nhmmer(taxa):
        _detection_legend(ax, "lower right", cfg)

    _draw_title(ax, gene_label, cfg)
    if show_meth:
        ax.text(0.02, 0.02, "NeighborNet outline (Bryant & Moulton 2004)",
                transform=ax.transAxes, fontsize=8, va="bottom", ha="left", color="#888")

    all_xe = np.array([x for x0,y0,x1,y1 in edges for x in (x0, x1)])
    all_ye = np.array([y for x0,y0,x1,y1 in edges for y in (y0, y1)])
    xpad = (all_xe.max() - all_xe.min()) * 0.20
    ypad = (all_ye.max() - all_ye.min()) * 0.20
    ax.set_xlim(all_xe.min() - xpad, all_xe.max() + xpad)
    ax.set_ylim(all_ye.min() - ypad, all_ye.max() + ypad)

    plt.tight_layout(pad=0.3)
    _save(fig, out_dir, f"{gene_label}_outline{suffix}")
    plt.close(fig)


# ─── Style 2: NJ + species bounding boxes ─────────────────────────────────────

def draw_nj_boxes(gene_label, gene, mldist_path, out_dir,
                  lookup, species2hex, chrom_lookup, detect_methods,
                  cfg=None, suffix=""):
    if not mldist_path.exists():
        print(f"  [SKIP] {mldist_path.name}")
        return

    taxa, D = parse_mldist(mldist_path)
    print(f"  [nj_boxes/{gene_label}] {len(taxa)} taxa")
    pos, edges = nj_layout(taxa, D)

    sp_taxa = defaultdict(list)
    for t in taxa:
        sp_taxa[lookup(t)].append(t)

    figsize   = tuple(_p(cfg, "figure", "figsize", default=(20, 20)))
    edge_col  = _p(cfg, "edges", "color",     default="#555555")
    edge_lw   = _p(cfg, "edges", "linewidth", default=0.7)
    s_blast   = _p(cfg, "markers", "size_blast",       default=55)
    s_nhmmer  = _p(cfg, "markers", "size_nhmmer",      default=80)
    ec_nhmmer = _p(cfg, "markers", "edge_color_nhmmer",default="#cc0000")
    lw_blast  = _p(cfg, "markers", "linewidth_blast",  default=0.4)
    lw_nhmmer = _p(cfg, "markers", "linewidth_nhmmer", default=1.0)
    lbl_fs    = _p(cfg, "labels", "fontsize", default=7.0)
    highlights = _p(cfg, "highlights", default={}) if cfg else {}
    show_meth = _p(cfg, "method_text", "show", default=True)
    box_lbl_fs = max(lbl_fs * 0.85, 12)   # species box label scales with leaf labels

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color=edge_col, linewidth=edge_lw, zorder=1)

    all_x  = np.array([pos[t][0] for t in taxa if t in pos])
    all_y  = np.array([pos[t][1] for t in taxa if t in pos])
    cx_all = all_x.mean()
    cy_all = all_y.mean()
    xrange = all_x.max() - all_x.min() if len(all_x) else 1.0

    PAD = 0.07
    box_info = {}
    for sp, members in sp_taxa.items():
        if sp == "Unknown":
            continue
        sp_hex = species2hex.get(sp, "#888888")
        pts_x  = [pos[t][0] for t in members if t in pos]
        pts_y  = [pos[t][1] for t in members if t in pos]
        if not pts_x:
            continue
        bx0 = min(pts_x) - PAD
        by0 = min(pts_y) - PAD * 0.8
        bx1 = max(pts_x) + PAD
        by1 = max(pts_y) + PAD * 0.8
        rect = plt.Rectangle((bx0, by0), bx1 - bx0, by1 - by0,
                              fill=True,
                              facecolor=(*hex_rgb(sp_hex), 0.10),
                              edgecolor=sp_hex, linewidth=2.0, zorder=2)
        ax.add_patch(rect)
        box_info[sp] = (bx0, by0, bx1, by1)

    for t in taxa:
        if t not in pos:
            continue
        p  = pos[t]
        sp = lookup(t)
        col = species2hex.get(sp, "#888888")
        label, strain, is_nhmmer = build_label(t, gene, chrom_lookup, detect_methods)

        hl       = highlights.get(strain, {})
        marker_s = hl.get("marker_size", s_nhmmer if is_nhmmer else s_blast)
        lbl_size = hl.get("label_fontsize", lbl_fs)
        lbl_fw   = hl.get("label_fontweight", "bold" if is_nhmmer else "normal")

        ax.scatter(p[0], p[1], c=col, s=marker_s,
                   edgecolors=ec_nhmmer if is_nhmmer else "white",
                   linewidths=lw_nhmmer if is_nhmmer else lw_blast,
                   zorder=4, marker="o")

        lx, ly, ha, va, angle = radial_label_pos(p[0], p[1], cx_all, cy_all, xrange)
        ax.text(lx, ly, label, ha=ha, va=va,
                fontsize=lbl_size, fontweight=lbl_fw,
                fontstyle="italic" if is_nhmmer else "normal",
                color="#cc0000" if is_nhmmer else "#111111",
                rotation=angle, rotation_mode="anchor",
                zorder=5)

    for sp, (bx0, by0, bx1, by1) in box_info.items():
        sp_hex = species2hex.get(sp, "#888888")
        parts  = sp.split()
        sp_lbl = f"S. {parts[1]}" if len(parts) >= 2 else sp
        box_cx = (bx0 + bx1) / 2
        box_cy = (by0 + by1) / 2
        dx = box_cx - cx_all
        dy = box_cy - cy_all
        if abs(dx) >= abs(dy):
            lx, ly, ha = (bx1 + 0.04, box_cy, "left") if dx >= 0 else (bx0 - 0.04, box_cy, "right")
        else:
            lx, ly, ha = (box_cx, by1 + 0.04, "center") if dy >= 0 else (box_cx, by0 - 0.04, "center")
        ax.text(lx, ly, sp_lbl, ha=ha, va="center",
                fontsize=box_lbl_fs, fontweight="bold", fontstyle="italic",
                color=sp_hex, zorder=6)

    if "extended" in gene_label.lower() and _has_nhmmer(taxa):
        _detection_legend(ax, "lower right", cfg)

    _draw_title(ax, gene_label, cfg)
    if show_meth:
        ax.text(0.02, 0.02, "NJ unrooted tree",
                transform=ax.transAxes, fontsize=8, va="bottom", ha="left", color="#888")

    tip_r  = np.sqrt((all_x - cx_all) ** 2 + (all_y - cy_all) ** 2)
    margin = np.median(tip_r) * 0.60
    ax.set_xlim(all_x.min() - margin, all_x.max() + margin)
    ax.set_ylim(all_y.min() - margin, all_y.max() + margin)

    plt.tight_layout(pad=0.3)
    _save(fig, out_dir, f"{gene_label}_nj_boxes{suffix}")
    plt.close(fig)


# ─── Style 3: NJ + colored leaves ─────────────────────────────────────────────

def draw_colored(gene_label, gene, mldist_path, out_dir,
                 lookup, species2hex, chrom_lookup, detect_methods,
                 cfg=None, suffix=""):
    if not mldist_path.exists():
        print(f"  [SKIP] {mldist_path.name}")
        return

    taxa, D = parse_mldist(mldist_path)
    print(f"  [colored/{gene_label}] {len(taxa)} taxa")
    pos, edges = nj_layout(taxa, D)

    figsize   = tuple(_p(cfg, "figure", "figsize", default=(20, 20)))
    edge_col  = _p(cfg, "edges", "color",     default="#aaaaaa")
    edge_lw   = _p(cfg, "edges", "linewidth", default=0.6)
    s_blast   = _p(cfg, "markers", "size_blast",       default=55)
    s_nhmmer  = _p(cfg, "markers", "size_nhmmer",      default=80)
    ec_nhmmer = _p(cfg, "markers", "edge_color_nhmmer",default="#cc0000")
    lw_blast  = _p(cfg, "markers", "linewidth_blast",  default=0.4)
    lw_nhmmer = _p(cfg, "markers", "linewidth_nhmmer", default=1.0)
    lbl_fs    = _p(cfg, "labels", "fontsize",    default=7.0)
    lbl_col_nm = _p(cfg, "labels", "color_nhmmer",default="#cc0000")
    highlights = _p(cfg, "highlights", default={}) if cfg else {}
    show_meth = _p(cfg, "method_text", "show", default=True)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color=edge_col, linewidth=edge_lw, zorder=1)

    all_x  = np.array([pos[t][0] for t in taxa if t in pos])
    all_y  = np.array([pos[t][1] for t in taxa if t in pos])
    cx     = all_x.mean()
    cy     = all_y.mean()
    xrange = all_x.max() - all_x.min() if len(all_x) else 1.0

    species_seen = {}
    for t in taxa:
        if t not in pos:
            continue
        p  = pos[t]
        sp = lookup(t)
        col = species2hex.get(sp, "#888888")
        label, strain, is_nhmmer = build_label(t, gene, chrom_lookup, detect_methods)

        hl       = highlights.get(strain, {})
        marker_s = hl.get("marker_size", s_nhmmer if is_nhmmer else s_blast)
        lbl_size = hl.get("label_fontsize", lbl_fs)
        lbl_fw   = hl.get("label_fontweight", "normal")

        ax.scatter(p[0], p[1], c=col, s=marker_s,
                   edgecolors=ec_nhmmer if is_nhmmer else "white",
                   linewidths=lw_nhmmer if is_nhmmer else lw_blast,
                   zorder=4, marker="o")

        lx, ly, ha, va, angle = radial_label_pos(p[0], p[1], cx, cy, xrange)
        ax.text(lx, ly, label, ha=ha, va=va,
                fontsize=lbl_size, fontweight=lbl_fw,
                fontstyle="italic" if is_nhmmer else "normal",
                color=lbl_col_nm if is_nhmmer else col,
                rotation=angle, rotation_mode="anchor",
                zorder=5)

        if sp not in species_seen:
            species_seen[sp] = col

    leg1 = _species_legend(ax, species_seen, "lower left", cfg)
    ax.add_artist(leg1)
    if "extended" in gene_label.lower() and _has_nhmmer(taxa):
        _detection_legend(ax, "lower right", cfg)

    _draw_title(ax, gene_label, cfg)
    if show_meth:
        ax.text(0.02, 0.02, "NJ unrooted tree",
                transform=ax.transAxes, fontsize=8, va="bottom", ha="left", color="#888")

    tip_r  = np.sqrt((all_x - cx) ** 2 + (all_y - cy) ** 2)
    margin = np.median(tip_r) * 0.60
    ax.set_xlim(all_x.min() - margin, all_x.max() + margin)
    ax.set_ylim(all_y.min() - margin, all_y.max() + margin)

    plt.tight_layout(pad=0.3)
    _save(fig, out_dir, f"{gene_label}_colored{suffix}")
    plt.close(fig)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    out_outline  = PROJECT / "resultados" / "Figuras_paper_2"
    out_nj_boxes = RESULTS / "networks_nj"
    out_colored  = RESULTS / "networks_colored"
    for d in (out_outline, out_nj_boxes, out_colored):
        d.mkdir(parents=True, exist_ok=True)

    print("[17] Loading metadata...")
    lookup, species2hex = load_metadata()
    detect_methods      = load_detect_methods()
    chrom_lookup        = load_chrom_lookup()
    cfg                 = load_style_config()

    for gene_label, (mldist_path, gene) in DATASETS.items():
        print(f"\n[17] === {gene_label} ===")
        draw_outline(gene_label, gene, mldist_path, out_outline,
                     lookup, species2hex, chrom_lookup, detect_methods,
                     cfg=cfg, suffix="_v2")
        draw_nj_boxes(gene_label, gene, mldist_path, out_nj_boxes,
                      lookup, species2hex, chrom_lookup, detect_methods,
                      cfg=cfg, suffix="_v2")
        draw_colored(gene_label, gene, mldist_path, out_colored,
                     lookup, species2hex, chrom_lookup, detect_methods,
                     cfg=cfg, suffix="_v2")

    print(f"\n[17] Done. Output: {RESULTS}/networks_{{outline,nj,colored}}/")

if __name__ == "__main__":
    main()
