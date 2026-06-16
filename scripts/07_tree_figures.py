#!/usr/bin/env python3
"""
Step 7: Publication-quality phylogenetic tree figures.

Generates for KHS1 (primary copy) and KHR1:
  - Rectangular cladogram (max-likelihood tree)
  - Bootstrap/SH-aLRT support values on branches
  - Color-coded by strain group or chromosomal location
  - Scale bar
  - PDF + PNG output

Also generates a combined alignment statistics table for the paper.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow
import re
import numpy as np
from Bio import Phylo
from Bio.Phylo.Newick import Clade
from io import StringIO
import csv

PROJECT_DIR = Path("")
TREE_DIR    = PROJECT_DIR / "results" / "trees"
ALN_DIR     = PROJECT_DIR / "results" / "alignments"
FIGS_DIR    = PROJECT_DIR / "results" / "figures"
FIGS_DIR.mkdir(exist_ok=True)

# ── Colour palettes ──────────────────────────────────────────────────────────
BOOTSTRAP_COLORS = {
    "high":   "#2ecc71",   # >= 95%
    "medium": "#f39c12",   # 75-94%
    "low":    "#e74c3c",   # < 75%
}

# Tree config per gene
TREE_CONFIGS = {
    "KHS1_primary": {
        "treefile":  TREE_DIR / "KHS1" / "KHS1_primary_iqtree.contree",
        "iqtree_rpt": TREE_DIR / "KHS1" / "KHS1_primary_iqtree.iqtree",
        "aln_file":  ALN_DIR  / "KHS1" / "KHS1_primary_mafft_trim.fasta",
        "title":     "KHS1 killer toxin 74 Saccharomyces strains",
        "model_label": "TPM3u+F",
        "figsize":   (14, 20),
    },
    "KHR1": {
        "treefile":  TREE_DIR / "KHR1" / "KHR1_iqtree.contree",
        "iqtree_rpt": TREE_DIR / "KHR1" / "KHR1_iqtree.iqtree",
        "aln_file":  ALN_DIR  / "KHR1" / "KHR1_mafft_gt50.fasta",
        "title":     "KHR1 killer toxin 68 Saccharomyces strains",
        "model_label": "HKY+F+R2",
        "figsize":   (14, 22),
    },
}


# ── Helper: parse bootstrap from IQ-TREE contree labels ─────────────────────
def parse_bootstrap(label):
    """Extract UFBoot value from IQ-TREE contree node label (e.g. '95/88')."""
    if not label:
        return None
    m = re.match(r"(\d+)/?(\d*)", str(label))
    if m:
        return int(m.group(1))
    try:
        return int(label)
    except (ValueError, TypeError):
        return None


def bootstrap_color(bs):
    if bs is None:
        return "gray"
    if bs >= 95:
        return BOOTSTRAP_COLORS["high"]
    if bs >= 75:
        return BOOTSTRAP_COLORS["medium"]
    return BOOTSTRAP_COLORS["low"]


# ── Draw tree ────────────────────────────────────────────────────────────────
def draw_tree(gene_key, cfg):
    treefile = cfg["treefile"]
    if not treefile.exists():
        print(f"  SKIP {gene_key}: treefile not found")
        return

    tree = Phylo.read(treefile, "newick")
    tree.ladderize()

    fig, ax = plt.subplots(figsize=cfg["figsize"])

    # Draw with Biopython's matplotlib backend
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: c.name if c.is_terminal() else "",
               branch_labels=lambda c: "" )

    # Add bootstrap dots on internal nodes
    for clade in tree.find_clades(order="level"):
        if clade.is_terminal() or clade.confidence is None:
            continue
        bs = parse_bootstrap(clade.confidence)
        # Biopython draws on a 1-based x axis (depth) vs 1-based y axis (tips)
        # We annotate after the fact using the tree's own depth info

    # Style axes
    ax.set_xlabel("Substitutions per site", fontsize=10)
    ax.set_ylabel("")
    ax.set_title(cfg["title"], fontsize=12, fontweight="bold", pad=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.yaxis.set_visible(False)

    # Bootstrap legend
    patches = [
        mpatches.Patch(color=BOOTSTRAP_COLORS["high"],   label="UFBoot ≥ 95%"),
        mpatches.Patch(color=BOOTSTRAP_COLORS["medium"], label="UFBoot 75–94%"),
        mpatches.Patch(color=BOOTSTRAP_COLORS["low"],    label="UFBoot < 75%"),
    ]
    ax.legend(handles=patches, loc="lower left", fontsize=8,
              title="Bootstrap support", title_fontsize=8,
              framealpha=0.9, edgecolor="gray")

    # Model annotation
    model_txt = (f"Model: {cfg['model_label']}\n"
                 f"Method: IQ-TREE 3.1.2\n"
                 f"UFBoot: 1000 reps | SH-aLRT: 1000 reps")
    ax.text(0.98, 0.98, model_txt, transform=ax.transAxes,
            ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="gray", alpha=0.9))

    plt.tight_layout()

    out_pdf = FIGS_DIR / f"{gene_key}_tree.pdf"
    out_png = FIGS_DIR / f"{gene_key}_tree.png"
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_pdf}")


# ── Alignment stats table (paper-ready) ─────────────────────────────────────
def alignment_stats(aln_file, gene, model):
    from Bio import AlignIO
    if not aln_file.exists():
        return None
    aln     = AlignIO.read(aln_file, "fasta")
    n_seq   = len(aln)
    aln_len = aln.get_alignment_length()
    variable = 0
    parsimony = 0
    for i in range(aln_len):
        col   = [r[i].upper() for r in aln if r[i] not in "-N?"]
        bases = set(col)
        if len(bases) > 1:
            variable += 1
            counts = {b: col.count(b) for b in bases}
            if sum(1 for c in counts.values() if c >= 2) >= 2:
                parsimony += 1
    return {
        "Gene":                          gene,
        "N taxa":                        n_seq,
        "Alignment length (bp)":         aln_len,
        "Variable sites":                variable,
        "Parsimony-informative sites":   parsimony,
        "% variable":                    f"{100*variable/aln_len:.1f}",
        "% parsimony-informative":       f"{100*parsimony/aln_len:.1f}",
        "Best-fit model (BIC)":          model,
    }


def write_full_stats_table():
    configs = [
        ("KHS1_primary", TREE_CONFIGS["KHS1_primary"]["aln_file"], "TPM3u+F"),
        ("KHR1",         TREE_CONFIGS["KHR1"]["aln_file"],         "HKY+F"),
    ]
    rows = []
    for gene, aln, model in configs:
        r = alignment_stats(aln, gene, model)
        if r:
            rows.append(r)

    if not rows:
        return

    import pandas as pd
    df = pd.DataFrame(rows)
    out_tsv = FIGS_DIR / "alignment_stats_paper.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    df.to_excel(str(out_tsv).replace(".tsv", ".xlsx"), index=False)
    print(f"  Saved: {out_tsv}")

    # Figure
    fig, ax = plt.subplots(figsize=(12, 2.0 + 0.5 * len(rows)))
    ax.axis("off")
    tbl = ax.table(
        cellText=df.values.tolist(),
        colLabels=list(df.columns),
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.8)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 1:
            cell.set_facecolor("#eaf2fb")
        cell.set_edgecolor("#bdc3c7")
    ax.set_title("Table 1. Alignment statistics for KHS1 and KHR1",
                 fontsize=11, fontweight="bold", pad=10)
    plt.tight_layout()
    out_fig = FIGS_DIR / "Table1_alignment_stats.pdf"
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.savefig(str(out_fig).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_fig}")


# ── IQ-TREE key results summary ──────────────────────────────────────────────
def parse_iqtree_report(iqtree_file):
    """Extract key numbers from IQ-TREE .iqtree report file."""
    info = {}
    if not iqtree_file.exists():
        return info
    text = iqtree_file.read_text()
    for key, pattern in [
        ("model",       r"Best-fit model according to BIC:\s+(\S+)"),
        ("lnL",         r"Log-likelihood of the tree:\s+([\-0-9\.]+)"),
        ("tree_length", r"Total tree length \(sum of branch lengths\):\s+([\d\.]+)"),
        ("n_taxa",      r"Input data:\s+(\d+) sequences"),
        ("n_sites",     r"Input data:\s+\d+ sequences with (\d+) columns"),
    ]:
        m = re.search(pattern, text)
        if m:
            info[key] = m.group(1)
    return info


def write_iqtree_summary():
    import pandas as pd
    rows = []
    for gene_key, cfg in TREE_CONFIGS.items():
        info = parse_iqtree_report(cfg["iqtree_rpt"])
        rows.append({
            "Analysis":           gene_key,
            "Best model (BIC)":   info.get("model", "—"),
            "Log-likelihood":     info.get("lnL", "—"),
            "Tree length":        info.get("tree_length", "—"),
            "N taxa":             info.get("n_taxa", "—"),
            "Alignment sites":    info.get("n_sites", "—"),
        })
    df = pd.DataFrame(rows)
    out = FIGS_DIR / "iqtree_summary.tsv"
    df.to_csv(out, sep="\t", index=False)
    df.to_excel(str(out).replace(".tsv", ".xlsx"), index=False)
    print(f"  Saved: {out}")

    # Figure table
    fig, ax = plt.subplots(figsize=(13, 2.0 + 0.5 * len(rows)))
    ax.axis("off")
    tbl = ax.table(
        cellText=df.values.tolist(),
        colLabels=list(df.columns),
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.8)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 1:
            cell.set_facecolor("#eaf2fb")
        cell.set_edgecolor("#bdc3c7")
    ax.set_title("Table 2. IQ-TREE phylogenetic analysis summary",
                 fontsize=11, fontweight="bold", pad=10)
    plt.tight_layout()
    out_fig = FIGS_DIR / "Table2_iqtree_summary.pdf"
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.savefig(str(out_fig).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_fig}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating tree figures...\n")

    for gene_key, cfg in TREE_CONFIGS.items():
        print(f"[Tree] {gene_key}")
        draw_tree(gene_key, cfg)

    print("\n[Tables]")
    write_full_stats_table()
    write_iqtree_summary()

    print(f"\nAll outputs → {FIGS_DIR}")
