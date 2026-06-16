#!/usr/bin/env python3

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
import matplotlib.lines as mlines
import dendropy
from dendropy.calculate import treecompare

PROJECT  = Path("")
RESULTS  = PROJECT / ""
OUT_DIR  = PROJECT / "" / ""
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENOMIC_TREE = RESULTS / "" / "supermatrix" / "phylogenomics_iqtree.contree"
STYLE_CONFIG_PATH = Path(__file__).parent / "style_config.json"

META_TSV = RESULTS / "01_metadata" / "genome_metadata.tsv"
SPP_CSV  = PROJECT / "Spp_Color_codes.csv"

GENE_TREES = {
    "KHS1":          RESULTS / "" / "KHS1" / "KHS1_primary_iqtree.contree",
    "KHR1":          RESULTS / "" / "KHR1" / "KHR1_iqtree.contree",
    "KHR1_extended": RESULTS / "05_hmm_discovery" / "candidatos_KHR1" / "KHR1_extended_iqtree.contree",
}

STRAIN_CORRECTIONS = {"UFRJ50816T": "UFRJ50816", "EM14S01-3B": "EM14S013B"}


# --- Data loaders -------------------------------------------------------------

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
        bn = name.split("__")[0] if "__" in name else name
        if bn in strain2species:
            return strain2species[bn]
        n = norm_strain(bn)
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


# --- Tree utilities -----------------------------------------------------------

def load_tree(path, taxon_namespace=None):
    kwargs = dict(path=str(path), schema="newick", preserve_underscores=True)
    if taxon_namespace is not None:
        kwargs["taxon_namespace"] = taxon_namespace
    return dendropy.Tree.get(**kwargs)

def leaf_labels(tree):
    return {nd.taxon.label for nd in tree.leaf_node_iter()}

def strip_nhmmer(name):
    return name.split("__")[0] if "__" in name else name

def collapse_extended_labels(tree):
    """Rename 'STRAIN__chrIX_..._nhmmer' -> 'STRAIN'.
    When a strain appears in both BLAST (plain) and nhmmer (extended) forms,
    keep the nhmmer copy (the novel sequence) and discard the BLAST copy.
    """
    nhmmer_strains = {nd.taxon.label.split("__")[0]
                      for nd in tree.leaf_node_iter()
                      if "__" in nd.taxon.label}
    to_prune = []
    for nd in tree.leaf_node_iter():
        if "__" in nd.taxon.label:
            nd.taxon.label = strip_nhmmer(nd.taxon.label)
        else:
            if nd.taxon.label in nhmmer_strains:
                to_prune.append(nd.taxon)
    for t in to_prune:
        tree.prune_taxa([t])
    tree.purge_taxon_namespace()

def prune_to_common(tree, keep_labels):
    to_prune = [nd.taxon for nd in tree.leaf_node_iter()
                if nd.taxon.label not in keep_labels]
    for taxon in to_prune:
        tree.prune_taxa([taxon])
    tree.purge_taxon_namespace()


# --- Robinson-Foulds ----------------------------------------------------------

def compute_rf(gene_label, gene_tree_path):
    if not gene_tree_path.exists():
        print(f"  [SKIP RF] {gene_tree_path.name} not found")
        return None

    # Load independently so label renaming doesn't corrupt a shared namespace
    g_tree = load_tree(GENOMIC_TREE)
    q_tree = load_tree(gene_tree_path)

    if "extended" in gene_label.lower():
        collapse_extended_labels(q_tree)

    common = leaf_labels(g_tree) & leaf_labels(q_tree)
    if len(common) < 4:
        print(f"  [SKIP RF] Too few common taxa: {len(common)}")
        return None

    prune_to_common(g_tree, common)
    prune_to_common(q_tree, common)

    # Merge into a shared namespace (required by treecompare)
    tns = dendropy.TaxonNamespace()
    g_tree.migrate_taxon_namespace(tns)
    q_tree.migrate_taxon_namespace(tns)

    g_tree.deroot()
    q_tree.deroot()
    g_tree.update_bipartitions()
    q_tree.update_bipartitions()

    rf    = treecompare.symmetric_difference(g_tree, q_tree)
    n     = len(common)
    max_rf = 2 * (n - 3) if n > 3 else 1
    rf_norm = rf / max_rf if max_rf > 0 else 0.0

    print(f"  RF {gene_label}: raw={rf}, n={n}, RF_norm={rf_norm:.3f} "
          f"({100*(1-rf_norm):.1f}% concordance)")
    return {"gene": gene_label, "n_common": n,
            "RF_raw": rf, "RF_max": max_rf, "RF_norm": round(rf_norm, 4)}


# --- Rectangular layout -------------------------------------------------------

def rectangular_layout(tree):
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


# --- Tanglegram ---------------------------------------------------------------

def draw_tanglegram(gene_label, gene_tree_path, lookup, species2hex, cfg):
    if not gene_tree_path.exists():
        print(f"  [SKIP tanglegram] {gene_tree_path.name} not found")
        return

    is_extended = "extended" in gene_label.lower()
    nhmmer_col  = _p(cfg, "labels", "color_nhmmer", default="#cc0000")

    g_tree = load_tree(GENOMIC_TREE)
    q_tree = load_tree(gene_tree_path)

    if is_extended:
        # Keep BOTH copies (BLAST plain + nhmmer extended) as separate leaves.
        # Prune genomic tree to strains that have at least one gene copy.
        common_strains = {strip_nhmmer(nd.taxon.label)
                          for nd in q_tree.leaf_node_iter()} & leaf_labels(g_tree)
        if len(common_strains) < 4:
            print(f"  [SKIP tanglegram] Too few common taxa: {len(common_strains)}")
            return
        prune_to_common(g_tree, common_strains)
        # Prune gene tree to taxa whose base strain is in common_strains
        to_prune = [nd.taxon for nd in q_tree.leaf_node_iter()
                    if strip_nhmmer(nd.taxon.label) not in common_strains]
        for t in to_prune:
            q_tree.prune_taxa([t])
        q_tree.purge_taxon_namespace()
        common = sorted(common_strains)
    else:
        common = sorted(leaf_labels(g_tree) & leaf_labels(q_tree))
        if len(common) < 4:
            print(f"  [SKIP tanglegram] Too few common taxa: {len(common)}")
            return
        prune_to_common(g_tree, set(common))
        prune_to_common(q_tree, set(common))

    g_tree.reroot_at_midpoint(update_bipartitions=True)
    q_tree.reroot_at_midpoint(update_bipartitions=True)

    gx, gy = rectangular_layout(g_tree)
    qx, qy = rectangular_layout(q_tree)
    gx_max = max(gx.values()) or 1.0
    qx_max = max(qx.values()) or 1.0

    edge_col = _p(cfg, "edges", "color",        default="#111111")
    edge_lw  = _p(cfg, "edges", "linewidth",    default=1.5)
    s_tip    = _p(cfg, "markers", "size_blast", default=120)
    lbl_fs   = _p(cfg, "labels", "fontsize",    default=14.0)

    g_leaf_y = {nd.taxon.label: gy[nd] for nd in g_tree.leaf_node_iter()}
    q_leaf_y = {nd.taxon.label: qy[nd] for nd in q_tree.leaf_node_iter()}
    gn = len(g_leaf_y)
    qn = len(q_leaf_y)

    n_rows = max(gn, qn)
    fig_h = max(20, n_rows * 0.32)
    fig, ax = plt.subplots(figsize=(30, fig_h))
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # X zones (axes coords): left tree, gap, right tree
    LX0, LX1 = 0.02, 0.38
    RX0, RX1 = 0.62, 0.98

    def mlx(raw):  return LX0 + (raw / gx_max) * (LX1 - LX0)
    def mrx(raw):  return RX1 - (raw / qx_max) * (RX1 - RX0)

    def ngy(y): return y / max(gn - 1, 1)
    def nqy(y): return y / max(qn - 1, 1)

    # Draw left tree (genomic)
    for nd in g_tree.preorder_node_iter():
        if nd.parent_node is None:
            continue
        ax.plot([mlx(gx[nd.parent_node]), mlx(gx[nd])],
                [ngy(gy[nd]), ngy(gy[nd])],
                color=edge_col, linewidth=edge_lw, transform=ax.transAxes, zorder=1)
    for nd in g_tree.preorder_node_iter():
        if nd.is_leaf():
            continue
        ch  = nd.child_nodes()
        ys  = [ngy(gy[c]) for c in ch]
        ax.plot([mlx(gx[nd])] * 2, [min(ys), max(ys)],
                color=edge_col, linewidth=edge_lw, transform=ax.transAxes, zorder=1)

    # Draw right tree (gene, mirrored)
    for nd in q_tree.preorder_node_iter():
        if nd.parent_node is None:
            continue
        ax.plot([mrx(qx[nd.parent_node]), mrx(qx[nd])],
                [nqy(qy[nd]), nqy(qy[nd])],
                color=edge_col, linewidth=edge_lw, transform=ax.transAxes, zorder=1)
    for nd in q_tree.preorder_node_iter():
        if nd.is_leaf():
            continue
        ch  = nd.child_nodes()
        ys  = [nqy(qy[c]) for c in ch]
        ax.plot([mrx(qx[nd])] * 2, [min(ys), max(ys)],
                color=edge_col, linewidth=edge_lw, transform=ax.transAxes, zorder=1)

    # Rank maps for concordance colouring
    g_sorted = sorted(g_leaf_y, key=lambda l: g_leaf_y[l])
    q_sorted = sorted(q_leaf_y, key=lambda l: q_leaf_y[l])
    g_rank = {lbl: i for i, lbl in enumerate(g_sorted)}
    q_rank = {lbl: i for i, lbl in enumerate(q_sorted)}

    # For extended: map strain → list of gene-tree labels (BLAST and/or nhmmer)
    if is_extended:
        strain_to_q: dict[str, list[str]] = {}
        for qlabel in q_leaf_y:
            strain_to_q.setdefault(strip_nhmmer(qlabel), []).append(qlabel)
    else:
        strain_to_q = {lbl: [lbl] for lbl in common}

    s = s_tip * 0.65
    species_seen: dict[str, str] = {}
    has_nhmmer = False

    for strain in common:
        if strain not in g_leaf_y:
            continue
        yg  = ngy(g_leaf_y[strain])
        sp  = lookup(strain)
        col = species2hex.get(sp, "#888888")

        # Left tip (genomic) — drawn once per strain
        ax.scatter(LX1, yg, c=col, s=s, edgecolors="white", linewidths=0.3,
                   transform=ax.transAxes, zorder=4, marker="o")
        ax.text(LX1 + 0.004, yg, strain, ha="left", va="center",
                fontsize=lbl_fs * 0.72, color=col,
                transform=ax.transAxes, zorder=5)

        if sp not in species_seen:
            species_seen[sp] = col

        # One or two gene-tree copies per strain
        for qlabel in strain_to_q.get(strain, []):
            if qlabel not in q_leaf_y:
                continue
            yq      = nqy(q_leaf_y[qlabel])
            is_nhmm = "__" in qlabel   # True only for nhmmer extended labels

            # Connecting line: concordant/discordant + solid/dashed for nhmmer
            conc     = abs(g_rank.get(strain, 0) - q_rank.get(qlabel, 0)) <= 2
            line_col = "#2a9d2a" if conc else "#cccccc"
            ax.plot([LX1, RX0], [yg, yq],
                    color=line_col, lw=0.9, alpha=0.75,
                    linestyle="--" if is_nhmm else "-",
                    transform=ax.transAxes, zorder=2)

            # Right tip
            if is_nhmm:
                tip_col    = nhmmer_col
                tip_marker = "^"
                tip_label  = f"{strain}*"
                has_nhmmer = True
            else:
                tip_col    = col
                tip_marker = "o"
                tip_label  = strain

            ax.scatter(RX0, yq, c=tip_col, s=s,
                       edgecolors=nhmmer_col if is_nhmm else "white",
                       linewidths=1.0 if is_nhmm else 0.3,
                       transform=ax.transAxes, zorder=4, marker=tip_marker)
            ax.text(RX0 - 0.004, yq, tip_label, ha="right", va="center",
                    fontsize=lbl_fs * 0.72, color=tip_col,
                    transform=ax.transAxes, zorder=5)

    # Titles
    fs = _p(cfg, "title", "fontsize", default=22)
    fw = _p(cfg, "title", "fontweight", default="bold")
    ax.text(LX0, 1.015, "Whole-genome phylogeny",
            transform=ax.transAxes, fontsize=fs * 0.45, fontweight=fw,
            va="bottom", ha="left")
    ax.text(RX1, 1.015, gene_label.replace("_extended", " (+nhmmer)"),
            transform=ax.transAxes, fontsize=fs * 0.45, fontweight=fw,
            fontstyle="italic", va="bottom", ha="right")

    # Legend
    leg_fs  = _p(cfg, "legend", "fontsize",       default=9)
    leg_tfs = _p(cfg, "legend", "title_fontsize", default=9)
    leg_nc  = _p(cfg, "legend", "ncol",           default=2)
    spp_patches = [
        mpatches.Patch(color=col,
                       label=f"S. {sp.split()[1]}" if len(sp.split()) >= 2 else sp)
        for sp, col in sorted(species_seen.items()) if sp != "Unknown"
    ]
    conc_lines = [
        mlines.Line2D([], [], color="#2a9d2a", lw=2, label="Concordant (rank ≤2)"),
        mlines.Line2D([], [], color="#cccccc", lw=2, label="Discordant"),
    ]
    if has_nhmmer:
        conc_lines += [
            mlines.Line2D([], [], color="#cccccc", lw=2, linestyle="--",
                          label="Discordant (2nd nhmmer locus)"),
            mlines.Line2D([], [], color=nhmmer_col, marker="^",
                          linestyle="None", markersize=8,
                          label="2nd KHR1 locus (nhmmer, ▲)"),
        ]
    leg_spp = ax.legend(handles=spp_patches, loc="lower left",
                        framealpha=0.92, edgecolor="#cccccc", ncol=leg_nc,
                        title="Species", title_fontsize=leg_tfs,
                        prop={"size": leg_fs, "style": "italic"},
                        bbox_to_anchor=(0.01, -0.015))
    ax.add_artist(leg_spp)
    ax.legend(handles=conc_lines, loc="lower right",
              fontsize=leg_fs, framealpha=0.92, edgecolor="#cccccc",
              title="Concordance", title_fontsize=leg_tfs,
              bbox_to_anchor=(0.99, -0.015))

    plt.tight_layout(pad=0.5)
    stem = f"tanglegram_{gene_label}_v2"
    for ext in ("pdf", "png", "svg"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else 150,
                    bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.close(fig)


# --- Main ---------------------------------------------------------------------

def main():
    print("[15] Loading metadata & style...")
    lookup, species2hex = load_metadata()
    cfg = load_style_config()

    print("\n[15] === Robinson-Foulds distances ===")
    rf_rows = []
    for gene_label, gene_tree_path in GENE_TREES.items():
        row = compute_rf(gene_label, gene_tree_path)
        if row:
            rf_rows.append(row)

    if rf_rows:
        df_rf = pd.DataFrame(rf_rows)
        out_rf = OUT_DIR / "RF_summary_v2.tsv"
        df_rf.to_csv(out_rf, sep="\t", index=False, float_format="%.4f")
        print(f"\n  RF summary:\n{df_rf.to_string(index=False)}")
        print(f"  Saved: {out_rf}")

    print("\n[15] === Tanglegrams ===")
    for gene_label, gene_tree_path in GENE_TREES.items():
        print(f"\n[15] {gene_label}")
        draw_tanglegram(gene_label, gene_tree_path, lookup, species2hex, cfg)

    print(f"\n[15] Done. Output: {OUT_DIR}")

if __name__ == "__main__":
    main()
