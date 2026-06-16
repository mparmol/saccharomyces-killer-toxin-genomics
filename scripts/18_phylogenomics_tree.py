#!/usr/bin/env python3
"""
18_phylogenomics_tree.py
Whole-genome ML phylogenomics tree — three layouts.

Input:  results_2/phylogenomics/supermatrix/phylogenomics_iqtree.contree
        (IQ-TREE UFBoot consensus tree; bootstrap values as node labels)
Output: results_2/phylogenomics/figures/{rectangular,circular,unrooted}_v2.{pdf,png,svg}
Style:  scripts/style_config.json
"""

import re
import json
import math
import numpy as np
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as mpatch
import dendropy
try:
    from scipy.spatial import ConvexHull as _ConvexHull
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

PROJECT  = Path("")
RESULTS  = PROJECT / ""
CONTREE  = RESULTS / "" / "supermatrix" / "phylogenomics_iqtree.contree"
OUT_DIR  = PROJECT / "results" / "paper_figures_2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

META_TSV = RESULTS / "01_metadata" / "genome_metadata.tsv"
SPP_CSV  = PROJECT / "Spp_Color_codes.csv"
STYLE_CONFIG_PATH = Path(__file__).parent / "style_config.json"

STRAIN_CORRECTIONS = {"UFRJ50816T": "UFRJ50816", "EM14S01-3B": "EM14S013B"}


# ─── Data loaders ─────────────────────────────────────────────────────────────

def norm_strain(s):
    return re.sub(r'[-_.\s]', '', s).upper()

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
            c = STRAIN_CORRECTIONS[s]
            strain2species[c] = sp
            norm2species[norm_strain(c)] = sp
    species2hex = {str(r["Species"]).strip(): str(r["HEX code"]).strip()
                   for _, r in spp.iterrows()}
    def lookup(name):
        if name in strain2species:
            return strain2species[name]
        n = norm_strain(name)
        if n in norm2species:
            return norm2species[n]
        if n and n[-1].isalpha() and n[:-1] in norm2species:
            return norm2species[n[:-1]]
        return "Unknown"
    return lookup, species2hex

def load_style_config():
    with open(STYLE_CONFIG_PATH) as f:
        return json.load(f)

def _p(cfg, *keys, default):
    if cfg is None:
        return default
    val = cfg
    for k in keys:
        if not isinstance(val, dict) or k not in val:
            return default
        val = val[k]
    return val

def parse_bootstrap(label):
    """Parse UFBoot value from node label; returns float or None."""
    if label is None:
        return None
    s = str(label).strip()
    try:
        # Handle 'X/Y' (alrt/ufboot) or plain 'X'
        return float(s.split('/')[-1])
    except ValueError:
        return None


# ─── Layout algorithms ────────────────────────────────────────────────────────

def rectangular_layout(tree):
    """
    Returns x_pos, y_pos dicts keyed by dendropy node.
    x = cumulative branch length from root.
    y = leaf index (0..n-1) in pre-order; internal = midpoint of children.
    """
    x_pos = {}
    def set_x(nd, depth):
        x_pos[nd] = depth
        for c in nd.child_nodes():
            set_x(c, depth + (c.edge_length or 0.0))
    set_x(tree.seed_node, 0.0)

    y_pos = {}
    counter = [0]
    def set_y(nd):
        if nd.is_leaf():
            y_pos[nd] = counter[0]
            counter[0] += 1
        else:
            for c in nd.child_nodes():
                set_y(c)
            ch = nd.child_nodes()
            y_pos[nd] = (y_pos[ch[0]] + y_pos[ch[-1]]) / 2
    set_y(tree.seed_node)
    return x_pos, y_pos


def circular_layout(tree):
    """
    Circular phylogram.
    angle = leaf index mapped to [0, 2π]; radius = cumulative branch length.
    """
    x_pos, y_pos_raw = rectangular_layout(tree)
    n_leaves = sum(1 for nd in tree.leaf_node_iter())

    angle_pos = {}
    for nd, y in y_pos_raw.items():
        angle_pos[nd] = y / max(n_leaves - 1, 1) * 2 * math.pi

    # Convert to Cartesian
    cart_pos = {}
    for nd in tree.preorder_node_iter():
        r = x_pos[nd]
        a = angle_pos[nd]
        cart_pos[nd] = np.array([r * math.cos(a), r * math.sin(a)])

    return x_pos, angle_pos, cart_pos


def unrooted_layout(tree):
    """Equal-angle unrooted layout."""
    lc = {}
    for nd in tree.postorder_node_iter():
        lc[nd] = 1 if nd.is_leaf() else sum(lc[c] for c in nd.child_nodes())

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

    pos = {tree.seed_node: np.zeros(2)}

    def place(nd, start, arc, pxy):
        mid = start + arc / 2
        el  = max(compress(nd.edge_length or 1e-6), 1e-4)
        xy  = np.array([pxy[0] + el * math.cos(mid),
                        pxy[1] + el * math.sin(mid)])
        pos[nd] = xy
        if nd.is_leaf():
            return
        children = list(nd.child_nodes())
        total = lc[nd]
        cur = start
        for ch in children:
            ch_arc = arc * lc[ch] / total
            place(ch, cur, ch_arc, xy)
            cur += ch_arc

    root = tree.seed_node
    children = list(root.child_nodes())
    total = lc[root]
    cur = 0.0
    for ch in children:
        place(ch, cur, 2 * math.pi * lc[ch] / total, (0.0, 0.0))
        cur += 2 * math.pi * lc[ch] / total

    tip_rs = [np.linalg.norm(pos[nd]) for nd in tree.leaf_node_iter() if nd in pos]
    scale  = np.percentile(tip_rs, 90) if tip_rs else 1.0
    scale  = scale if scale > 1e-9 else 1.0
    pos = {nd: p / scale for nd, p in pos.items()}

    edges = [(pos[nd.parent_node][0], pos[nd.parent_node][1],
              pos[nd][0], pos[nd][1])
             for nd in tree.preorder_node_iter()
             if nd.parent_node is not None and nd in pos and nd.parent_node in pos]
    return pos, edges


# ─── Shared drawing helpers ────────────────────────────────────────────────────

BOOT_HIGH  = 95   # UFBoot threshold — strong support (filled circle)
BOOT_MED   = 70   # UFBoot threshold — moderate support (open circle)

def _boot_markers(ax, nodes, x_pos, y_pos, boot_high=BOOT_HIGH, boot_med=BOOT_MED):
    """Draw bootstrap circles at internal nodes."""
    for nd in nodes:
        if nd.is_leaf():
            continue
        val = parse_bootstrap(nd.label)
        if val is None or val < boot_med:
            continue
        x, y = x_pos[nd], y_pos[nd]
        if val >= boot_high:
            ax.scatter(x, y, s=60, c="#1a7a1a", zorder=6, marker="o", linewidths=0)
        else:
            ax.scatter(x, y, s=40, c="none", edgecolors="#1a7a1a",
                       linewidths=1.2, zorder=6, marker="o")

def _boot_markers_xy(ax, nodes, pos_dict, boot_high=BOOT_HIGH, boot_med=BOOT_MED):
    """Bootstrap markers using a pos dict keyed by node → (x, y) array."""
    for nd in nodes:
        if nd.is_leaf() or nd not in pos_dict:
            continue
        val = parse_bootstrap(nd.label)
        if val is None or val < boot_med:
            continue
        x, y = pos_dict[nd]
        if val >= boot_high:
            ax.scatter(x, y, s=60, c="#1a7a1a", zorder=6, marker="o", linewidths=0)
        else:
            ax.scatter(x, y, s=40, c="none", edgecolors="#1a7a1a",
                       linewidths=1.2, zorder=6, marker="o")

def _species_legend(ax, species_seen, cfg, loc="lower left"):
    fs  = _p(cfg, "legend", "fontsize",       default=9)
    tfs = _p(cfg, "legend", "title_fontsize", default=9)
    nc  = _p(cfg, "legend", "ncol",           default=2)
    patches = [
        mpatches.Patch(color=col,
                       label=f"S. {sp.split()[1]}" if len(sp.split()) >= 2 else sp)
        for sp, col in sorted(species_seen.items()) if sp != "Unknown"
    ]
    boot_patches = [
        mpatches.Patch(facecolor="#1a7a1a", edgecolor="none",   label=f"UFBoot ≥{BOOT_HIGH}"),
        mpatches.Patch(facecolor="none",    edgecolor="#1a7a1a", label=f"UFBoot {BOOT_MED}–{BOOT_HIGH-1}"),
    ]
    leg1 = ax.legend(handles=patches, loc=loc, framealpha=0.92, edgecolor="#cccccc",
                     ncol=nc, title="Species", title_fontsize=tfs,
                     prop={"size": fs, "style": "italic"})
    ax.add_artist(leg1)
    ax.legend(handles=boot_patches, loc="lower right", fontsize=max(fs * 0.8, 8),
              framealpha=0.92, edgecolor="#cccccc",
              title="Bootstrap", title_fontsize=max(tfs * 0.8, 8))

def _draw_title(ax, title, cfg):
    fs = _p(cfg, "title", "fontsize",   default=22)
    fw = _p(cfg, "title", "fontweight", default="bold")
    ax.text(0.02, 0.98, title,
            transform=ax.transAxes, fontsize=fs, fontweight=fw,
            va="top", ha="left")

def _save(fig, stem):
    for ext in ("pdf", "png", "svg"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else 150,
                    bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.close(fig)


# ─── Species annotation helpers ──────────────────────────────────────────────

def _sp_short(sp):
    parts = sp.split()
    return f"S. {parts[1]}" if len(parts) >= 2 else sp


def _species_brackets_rect(ax, tree, x_pos, y_pos, lookup, species2hex, x_max, lbl_fs):
    """Vertical bracket lines to the right of the rectangular tree, one per species."""
    sp_ys = {}
    for nd in tree.leaf_node_iter():
        sp = lookup(nd.taxon.label)
        if sp == "Unknown":
            continue
        sp_ys.setdefault(sp, []).append(y_pos[nd])

    bx       = x_max * 1.24
    tick_len = x_max * 0.012
    lx       = bx + tick_len * 0.9

    for sp, ys in sp_ys.items():
        col  = species2hex.get(sp, "#888888")
        y0   = min(ys) - 0.35
        y1   = max(ys) + 0.35
        # Vertical bar
        ax.plot([bx, bx], [y0, y1],
                color=col, linewidth=3.5, solid_capstyle="butt", zorder=4)
        # Short horizontal end-ticks
        for y_end in (y0, y1):
            ax.plot([bx - tick_len, bx], [y_end, y_end],
                    color=col, linewidth=2.5, zorder=4)
        # Label
        ax.text(lx, (y0 + y1) / 2, _sp_short(sp),
                ha="left", va="center", fontsize=lbl_fs * 1.76,
                color=col, fontstyle="italic", fontweight="bold", zorder=5)


def _species_brackets_circ(ax, tree, x_pos, angle_pos, lookup, species2hex,
                            x_max, lbl_fs):
    """Arc brackets just outside the circular tree, one per species."""
    sp_angles = {}
    for nd in tree.leaf_node_iter():
        sp = lookup(nd.taxon.label)
        if sp == "Unknown":
            continue
        sp_angles.setdefault(sp, []).append(angle_pos[nd])

    br = x_max * 1.40   # bracket radius
    lr = x_max * 1.50   # label radius

    for sp, angles in sp_angles.items():
        col  = species2hex.get(sp, "#888888")
        a_s  = sorted(angles)
        a0, a1 = a_s[0], a_s[-1]

        # Handle wraparound: if the span > π the arc passes through 0
        if a1 - a0 > math.pi:
            a0, a1 = a1, a_s[0] + 2 * math.pi   # short arc going through 0

        # Inset the drawn arc slightly so adjacent-species endpoints don't touch
        gap = min(0.05, abs(a1 - a0) * 0.08)
        a0d, a1d = a0 + gap, a1 - gap   # drawn arc limits

        n_seg = max(4, int(abs(a1d - a0d) / 0.04))
        ts = np.linspace(a0d, a1d, n_seg)
        ax.plot(br * np.cos(ts), br * np.sin(ts),
                color=col, linewidth=3.5, solid_capstyle="round", zorder=4)

        # End ticks at the inset arc ends
        for a_end in (a0d, a1d):
            ax.plot([br * 0.975 * math.cos(a_end), br * 1.025 * math.cos(a_end)],
                    [br * 0.975 * math.sin(a_end), br * 1.025 * math.sin(a_end)],
                    color=col, linewidth=2.5, zorder=4)

        # Rotated label at arc midpoint
        a_mid = (a0d + a1d) / 2
        lx = lr * math.cos(a_mid)
        ly = lr * math.sin(a_mid)
        ha  = "left" if math.cos(a_mid) > 0.1 else ("right" if math.cos(a_mid) < -0.1 else "center")
        rot = math.degrees(a_mid)
        if math.cos(a_mid) < -0.1:
            rot += 180
        ax.text(lx, ly, _sp_short(sp), ha=ha, va="center",
                fontsize=lbl_fs * 1.76, color=col, fontstyle="italic", fontweight="bold",
                rotation=rot, rotation_mode="anchor", zorder=6)


def _species_hulls_unrooted(ax, sp_pts, all_cx, all_cy, species2hex, lbl_fs, xrange):
    """Convex-hull polygon per species around unrooted fan layout.

    sp_pts: dict {species: [np.array, ...]} — flat list of points that should be
    enclosed (both node positions AND estimated label-end positions).
    """
    pad = xrange * 0.040  # pad beyond hull vertices

    for sp, pts_list in sp_pts.items():
        col = species2hex.get(sp, "#888888")
        pts = np.array(pts_list)
        # Deduplicate (keeps unique rows; order unimportant for hull)
        pts = np.unique(pts.round(8), axis=0)

        if len(pts) <= 1:
            angles_h = np.linspace(0, 2 * math.pi, 7)[:-1]
            poly_v = pts[0] + pad * 2 * np.column_stack(
                [np.cos(angles_h), np.sin(angles_h)])
        elif len(pts) == 2:
            p1, p2 = pts[0], pts[1]
            d = p2 - p1
            dist = np.linalg.norm(d)
            u = d / dist if dist > 1e-9 else np.array([1.0, 0.0])
            perp = np.array([-u[1], u[0]])
            poly_v = np.array([
                p1 - u * pad + perp * pad,
                p2 + u * pad + perp * pad,
                p2 + u * pad - perp * pad,
                p1 - u * pad - perp * pad,
            ])
        else:
            if _HAS_SCIPY:
                try:
                    hull = _ConvexHull(pts)
                    hull_v = pts[hull.vertices]
                except Exception:
                    c = pts.mean(axis=0)
                    angs = np.arctan2(pts[:, 1] - c[1], pts[:, 0] - c[0])
                    hull_v = pts[np.argsort(angs)]
            else:
                c = pts.mean(axis=0)
                angs = np.arctan2(pts[:, 1] - c[1], pts[:, 0] - c[0])
                hull_v = pts[np.argsort(angs)]

            centroid = hull_v.mean(axis=0)
            dirs = hull_v - centroid
            norms = np.linalg.norm(dirs, axis=1, keepdims=True)
            norms = np.where(norms < 1e-9, 1.0, norms)
            poly_v = centroid + dirs + (dirs / norms) * pad

        polygon = mpatch.Polygon(poly_v, closed=True, fill=False,
                                 edgecolor=col, linewidth=2.5, zorder=3,
                                 linestyle="-")
        ax.add_patch(polygon)

        # Label: project away from tree centre
        centroid = poly_v.mean(axis=0)
        dx, dy = centroid[0] - all_cx, centroid[1] - all_cy
        dist = math.hypot(dx, dy)
        nx, ny = (dx / dist, dy / dist) if dist > 1e-9 else (1.0, 0.0)
        sp_r  = np.linalg.norm(poly_v - centroid, axis=1).max()
        lx = centroid[0] + nx * (sp_r + pad * 1.8)
        ly = centroid[1] + ny * (sp_r + pad * 1.8)
        ha = "left" if nx > 0.1 else ("right" if nx < -0.1 else "center")
        ax.text(lx, ly, _sp_short(sp), ha=ha, va="center",
                fontsize=lbl_fs * 1.76, color=col,
                fontstyle="italic", fontweight="bold", zorder=6)


# ─── Layout 1: Rectangular phylogram ─────────────────────────────────────────

def draw_rectangular(tree, lookup, species2hex, cfg):
    x_pos, y_pos = rectangular_layout(tree)
    n_leaves = sum(1 for _ in tree.leaf_node_iter())

    edge_col  = _p(cfg, "edges", "color",     default="#111111")
    edge_lw   = _p(cfg, "edges", "linewidth", default=1.5)
    s_tip     = _p(cfg, "markers", "size_blast", default=120)
    lbl_fs    = _p(cfg, "labels", "fontsize",    default=14.0)

    fig_w = 20
    fig_h = max(20, n_leaves * 0.32)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    fig.patch.set_facecolor("white")

    x_max = max(x_pos.values())

    # Draw branches
    for nd in tree.preorder_node_iter():
        if nd.parent_node is None:
            continue
        xn, yn = x_pos[nd], y_pos[nd]
        xp, yp = x_pos[nd.parent_node], y_pos[nd.parent_node]
        # Horizontal branch
        ax.plot([xp, xn], [yn, yn], color=edge_col, linewidth=edge_lw, zorder=1)

    # Vertical connector at each internal node
    for nd in tree.preorder_node_iter():
        if nd.is_leaf():
            continue
        children = nd.child_nodes()
        xn = x_pos[nd]
        ys = [y_pos[c] for c in children]
        ax.plot([xn, xn], [min(ys), max(ys)], color=edge_col, linewidth=edge_lw, zorder=1)

    # Bootstrap markers and tip dots + labels
    species_seen = {}
    for nd in tree.preorder_node_iter():
        xn, yn = x_pos[nd], y_pos[nd]
        if nd.is_leaf():
            sp  = lookup(nd.taxon.label)
            col = species2hex.get(sp, "#888888")
            ax.scatter(xn, yn, c=col, s=s_tip, zorder=4,
                       edgecolors="white", linewidths=0.4, marker="o")
            gap = x_max * 0.008
            ax.text(xn + gap, yn, nd.taxon.label,
                    ha="left", va="center", fontsize=lbl_fs,
                    color=col, zorder=5)
            if sp not in species_seen:
                species_seen[sp] = col
        else:
            val = parse_bootstrap(nd.label)
            if val is not None and val >= BOOT_MED:
                c = "#1a7a1a" if val >= BOOT_HIGH else "none"
                ec = "#1a7a1a"
                ax.scatter(xn, yn, s=60, c=c, edgecolors=ec,
                           linewidths=1.2, zorder=6, marker="o")

    # Scale bar (0.01 substitutions/site)
    bar_len = 0.01
    ax.plot([0, bar_len], [-1.5, -1.5], color=edge_col, linewidth=edge_lw * 1.5)
    ax.text(bar_len / 2, -2.2, "0.01 sub/site",
            ha="center", va="top", fontsize=lbl_fs * 0.85, color="#333")

    _species_brackets_rect(ax, tree, x_pos, y_pos, lookup, species2hex, x_max, lbl_fs)

    ax.set_xlim(-x_max * 0.02, x_max * 1.72)
    ax.set_ylim(-3, n_leaves + 1)

    _species_legend(ax, species_seen, cfg, "lower left")
    _draw_title(ax, "Phylogenomics", cfg)
    ax.text(0.02, 0.01, "ML tree (IQ-TREE, MFP+MERGE, UFBoot 1000)",
            transform=ax.transAxes, fontsize=max(lbl_fs * 0.75, 9),
            va="bottom", ha="left", color="#888")

    plt.tight_layout(pad=0.3)
    _save(fig, "phylogenomics_rectangular_v2")


# ─── Layout 2: Circular phylogram ────────────────────────────────────────────

def draw_circular(tree, lookup, species2hex, cfg):
    x_pos, angle_pos, cart_pos = circular_layout(tree)

    # ── Non-uniform spacing: add angular gap between species groups ────────────
    leaves_ord = [nd for nd in tree.preorder_node_iter() if nd.is_leaf()]
    sp_seq     = [lookup(nd.taxon.label) for nd in leaves_ord]
    n_lv       = len(leaves_ord)
    GAP        = 0.8   # each interspecies gap = 0.8 normal leaf slots

    n_trans = sum(sp_seq[i] != sp_seq[i - 1] for i in range(1, n_lv))
    if n_lv > 0 and sp_seq[0] != sp_seq[-1]:
        n_trans += 1   # wraparound gap

    slot = 2 * math.pi / max((n_lv - 1) + n_trans * GAP, 1)
    nl_a = {}   # new leaf angle
    cum  = 0.0
    for i, nd in enumerate(leaves_ord):
        nl_a[nd] = cum
        if i < n_lv - 1:
            cum += slot * (1 + GAP) if sp_seq[i + 1] != sp_seq[i] else slot

    # Propagate to internal nodes: mean angle of leaf descendants (bottom-up)
    node_lf_angles = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            node_lf_angles[nd] = [nl_a[nd]] if nd in nl_a else []
            if nd in nl_a:
                angle_pos[nd] = nl_a[nd]
        else:
            lf_a = []
            for c in nd.child_nodes():
                lf_a.extend(node_lf_angles.get(c, []))
            node_lf_angles[nd] = lf_a
            if lf_a:
                angle_pos[nd] = sum(lf_a) / len(lf_a)

    # Refresh Cartesian positions with new angles
    for nd in tree.preorder_node_iter():
        if nd in x_pos and nd in angle_pos:
            r = x_pos[nd]
            a = angle_pos[nd]
            cart_pos[nd] = np.array([r * math.cos(a), r * math.sin(a)])
    # ──────────────────────────────────────────────────────────────────────────

    edge_col  = _p(cfg, "edges", "color",     default="#111111")
    edge_lw   = _p(cfg, "edges", "linewidth", default=1.5)
    s_tip     = _p(cfg, "markers", "size_blast", default=120)
    lbl_fs    = _p(cfg, "labels", "fontsize",    default=14.0)

    fig, ax = plt.subplots(figsize=(26, 26))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    x_max = max(x_pos.values())
    n_seg = 60   # arc approximation segments

    # Draw branches
    for nd in tree.preorder_node_iter():
        if nd.parent_node is None:
            continue
        xn = x_pos[nd]
        xp = x_pos[nd.parent_node]
        an = angle_pos[nd]
        # Horizontal (radial) branch
        ax.plot([xp * math.cos(an), xn * math.cos(an)],
                [xp * math.sin(an), xn * math.sin(an)],
                color=edge_col, linewidth=edge_lw, zorder=1)

    # Vertical (arc) connectors at internal nodes
    for nd in tree.preorder_node_iter():
        if nd.is_leaf():
            continue
        xn = x_pos[nd]
        angles = sorted(angle_pos[c] for c in nd.child_nodes())
        a_min, a_max = angles[0], angles[-1]
        # Draw arc from a_min to a_max at radius xn
        if a_max - a_min > math.pi:   # take shorter arc
            a_min, a_max = a_max, a_min + 2 * math.pi
        ts = np.linspace(a_min, a_max, n_seg)
        ax.plot(xn * np.cos(ts), xn * np.sin(ts),
                color=edge_col, linewidth=edge_lw, zorder=1)

    species_seen = {}
    for nd in tree.preorder_node_iter():
        if nd not in cart_pos:
            continue
        xy = cart_pos[nd]
        if nd.is_leaf():
            sp  = lookup(nd.taxon.label)
            col = species2hex.get(sp, "#888888")
            ax.scatter(xy[0], xy[1], c=col, s=s_tip, zorder=4,
                       edgecolors="white", linewidths=0.4, marker="o")
            # Radial label
            an = angle_pos[nd]
            r  = x_max * 1.07
            lx, ly = r * math.cos(an), r * math.sin(an)
            ha = "left" if -math.pi / 2 <= an <= math.pi / 2 else "right"
            rot = math.degrees(an)
            if ha == "right":
                rot += 180
            ax.text(lx, ly, nd.taxon.label,
                    ha=ha, va="center", fontsize=lbl_fs,
                    rotation=rot, rotation_mode="anchor",
                    color=col, zorder=5)
            if sp not in species_seen:
                species_seen[sp] = col
        else:
            val = parse_bootstrap(nd.label)
            if val is not None and val >= BOOT_MED:
                c  = "#1a7a1a" if val >= BOOT_HIGH else "none"
                ec = "#1a7a1a"
                ax.scatter(xy[0], xy[1], s=60, c=c, edgecolors=ec,
                           linewidths=1.2, zorder=6, marker="o")

    _species_brackets_circ(ax, tree, x_pos, angle_pos, lookup, species2hex,
                            x_max, lbl_fs)

    pad = x_max * 2.05
    ax.set_xlim(-pad, pad)
    ax.set_ylim(-pad, pad)

    _species_legend(ax, species_seen, cfg, "lower left")
    _draw_title(ax, "Phylogenomics", cfg)
    plt.tight_layout(pad=0.3)
    _save(fig, "phylogenomics_circular_v2")


# ─── Layout 3: Unrooted fan ───────────────────────────────────────────────────

def draw_unrooted(tree_unrooted, lookup, species2hex, cfg):
    pos, edges = unrooted_layout(tree_unrooted)

    edge_col  = _p(cfg, "edges", "color",     default="#111111")
    edge_lw   = _p(cfg, "edges", "linewidth", default=1.5)
    s_tip     = _p(cfg, "markers", "size_blast", default=120)
    lbl_fs    = _p(cfg, "labels", "fontsize",    default=14.0)

    fig, ax = plt.subplots(figsize=(24, 24))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for x0, y0, x1, y1 in edges:
        ax.plot([x0, x1], [y0, y1], color=edge_col, linewidth=edge_lw, zorder=1)

    all_lx = np.array([pos[nd][0] for nd in tree_unrooted.leaf_node_iter() if nd in pos])
    all_ly = np.array([pos[nd][1] for nd in tree_unrooted.leaf_node_iter() if nd in pos])
    cx, cy = all_lx.mean(), all_ly.mean()
    xrange = all_lx.max() - all_lx.min() if len(all_lx) else 1.0

    off = xrange * 0.025

    sp_hull_pts = {}   # species -> list of node positions (dots only)
    species_seen = {}

    for nd in tree_unrooted.preorder_node_iter():
        if nd not in pos:
            continue
        p = pos[nd]
        if nd.is_leaf():
            sp  = lookup(nd.taxon.label)
            col = species2hex.get(sp, "#888888")
            ax.scatter(p[0], p[1], c=col, s=s_tip, zorder=4,
                       edgecolors="white", linewidths=0.4, marker="o")

            dx, dy = p[0] - cx, p[1] - cy
            dist   = math.hypot(dx, dy)
            nx, ny = (dx / dist, dy / dist) if dist > 1e-9 else (1.0, 0.0)
            lx, ly = p[0] + nx * off, p[1] + ny * off
            angle = math.degrees(math.atan2(dy, dx))
            if nx < -0.1:
                angle += 180
                ha = "right"
            elif nx > 0.1:
                ha = "left"
            else:
                ha = "center"
            ax.text(lx, ly, nd.taxon.label, ha=ha, va="center", fontsize=lbl_fs,
                    color=col, zorder=5, rotation=angle, rotation_mode="anchor")

            sp_hull_pts.setdefault(sp, []).append(p.copy())

            if sp not in species_seen:
                species_seen[sp] = col
        else:
            val = parse_bootstrap(nd.label)
            if val is not None and val >= BOOT_MED and nd in pos:
                c  = "#1a7a1a" if val >= BOOT_HIGH else "none"
                ax.scatter(p[0], p[1], s=60, c=c, edgecolors="#1a7a1a",
                           linewidths=1.2, zorder=6, marker="o")

    _species_hulls_unrooted(ax, sp_hull_pts, cx, cy, species2hex, lbl_fs, xrange)

    tip_r  = np.sqrt((all_lx - cx)**2 + (all_ly - cy)**2)
    margin = np.median(tip_r) * 1.50   # room for rotated labels + species name labels
    ax.set_xlim(all_lx.min() - margin, all_lx.max() + margin)
    ax.set_ylim(all_ly.min() - margin, all_ly.max() + margin)

    _species_legend(ax, species_seen, cfg, "lower left")
    _draw_title(ax, "Phylogenomics", cfg)
    plt.tight_layout(pad=0.3)
    _save(fig, "phylogenomics_unrooted_v2")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("[18] Loading tree...")
    tree = dendropy.Tree.get(
        path=str(CONTREE),
        schema="newick",
        preserve_underscores=True,
    )
    tree.reroot_at_midpoint(update_bipartitions=True)
    n = sum(1 for _ in tree.leaf_node_iter())
    print(f"[18] {n} taxa loaded")

    # Unrooted copy (before midpoint root alters topology for fan layout)
    tree_unrooted = dendropy.Tree.get(
        path=str(CONTREE),
        schema="newick",
        preserve_underscores=True,
    )

    print("[18] Loading metadata & style...")
    lookup, species2hex = load_metadata()
    cfg = load_style_config()

    print("\n[18] === Rectangular ===")
    draw_rectangular(tree, lookup, species2hex, cfg)

    print("\n[18] === Circular ===")
    draw_circular(tree, lookup, species2hex, cfg)

    print("\n[18] === Unrooted fan ===")
    draw_unrooted(tree_unrooted, lookup, species2hex, cfg)

    print(f"\n[18] Done. Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
