#!/usr/bin/env python3
"""
Step 8: Extra publication figures.
  1. Chromosomal location of KHS1 and KHR1 hits across strains (dot plot)
  2. Pairwise % identity heatmap (from IQ-TREE ML distances)
  3. KHR1 bootstrap support summary plot
"""

from pathlib import Path
import csv
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

PROJECT_DIR = Path("")
RESULTS_DIR = PROJECT_DIR / "results"
FIGS_DIR    = RESULTS_DIR / "figures"
FIGS_DIR.mkdir(exist_ok=True)
SEQ_DIR     = RESULTS_DIR / "sequences"
TREE_DIR    = RESULTS_DIR / "trees"

GENES = ["KHS1", "KHR1"]

# Roman to int for chromosome sorting
ROMAN = {"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8,
         "IX":9,"X":10,"XI":11,"XII":12,"XIII":13,"XIV":14,"XV":15,"XVI":16}

def roman_sort_key(chrom):
    """Return integer sort key for 'chrXIV' style chromosome names."""
    c = chrom.replace("chr","").upper()
    return ROMAN.get(c, 99)


# ════════════════════════════════════════════════════════════════════════════
# 1. Chromosomal location dot plot
# ════════════════════════════════════════════════════════════════════════════
def plot_chromosomal_location(gene):
    detail_file = SEQ_DIR / gene / f"{gene}_blast_hits_detail.tsv"
    if not detail_file.exists():
        return

    df = pd.read_csv(detail_file, sep="\t")
    # Only primary copy (copy==1)
    df = df[df["copy"] == 1].copy()

    # Sort strains alphabetically
    strains = sorted(df["strain"].unique())
    chroms  = sorted(df["chrom"].unique(), key=roman_sort_key)

    # Map to integer positions
    strain_idx = {s: i for i, s in enumerate(strains)}
    chrom_idx  = {c: i for i, c in enumerate(chroms)}

    fig, ax = plt.subplots(figsize=(max(6, len(chroms)*0.6),
                                    max(8, len(strains)*0.22)))

    # Colour by strand
    strand_colors = {"+": "#3498db", "-": "#e74c3c"}

    for _, row in df.iterrows():
        x = chrom_idx.get(row["chrom"], 0)
        y = strain_idx.get(row["strain"], 0)
        c = strand_colors.get(row["strand"], "gray")
        ax.scatter(x, y, color=c, s=60, zorder=3, alpha=0.85,
                   edgecolors="white", linewidths=0.4)

    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(strains)))
    ax.set_yticklabels(strains, fontsize=7)
    ax.set_xlim(-0.5, len(chroms) - 0.5)
    ax.set_ylim(-0.5, len(strains) - 0.5)
    ax.set_xlabel("Chromosome", fontsize=10)
    ax.set_ylabel("Strain", fontsize=10)
    ax.set_title(f"Chromosomal location of {gene}\n(primary copy per strain)",
                 fontsize=11, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.3, zorder=0)

    patches = [
        mpatches.Patch(color=strand_colors["+"], label="(+) strand"),
        mpatches.Patch(color=strand_colors["-"], label="(−) strand"),
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=8)

    plt.tight_layout()
    out = FIGS_DIR / f"{gene}_chromosomal_location.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════════
# 2. Pairwise ML distance heatmap (from IQ-TREE .mldist)
# ════════════════════════════════════════════════════════════════════════════
def plot_mldist_heatmap(gene, tree_prefix):
    mldist_file = TREE_DIR / gene / f"{tree_prefix}.mldist"
    if not mldist_file.exists():
        print(f"  SKIP mldist: {mldist_file} not found")
        return

    lines = mldist_file.read_text().strip().splitlines()
    n     = int(lines[0].strip())
    taxa  = []
    mat   = []
    for line in lines[1:n+1]:
        parts = line.split()
        taxa.append(parts[0])
        mat.append([float(x) for x in parts[1:]])

    df  = pd.DataFrame(mat, index=taxa, columns=taxa)

    # Limit to manageable size for display (cluster if > 60)
    if len(taxa) > 70:
        # Show only the first 60 taxa for clarity
        df = df.iloc[:60, :60]

    fig_size = max(10, len(df) * 0.2)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    cmap = LinearSegmentedColormap.from_list(
        "dist", ["#2ecc71", "#f7dc6f", "#e74c3c"])

    sns.heatmap(df, ax=ax, cmap=cmap, square=True,
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "ML distance", "shrink": 0.6},
                linewidths=0, annot=False)

    tick_fontsize = max(4, min(8, 200 // len(df)))
    ax.tick_params(axis="x", labelsize=tick_fontsize, rotation=90)
    ax.tick_params(axis="y", labelsize=tick_fontsize, rotation=0)
    ax.set_title(f"{gene} — pairwise ML distances",
                 fontsize=11, fontweight="bold", pad=10)

    plt.tight_layout()
    out = FIGS_DIR / f"{gene}_mldist_heatmap.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════════
# 3. % Identity barplot per strain (relative to query)
# ════════════════════════════════════════════════════════════════════════════
def plot_identity_barplot(gene):
    detail_file = SEQ_DIR / gene / f"{gene}_blast_hits_detail.tsv"
    if not detail_file.exists():
        return

    df = pd.read_csv(detail_file, sep="\t")
    df = df[df["copy"] == 1].copy()
    df = df.sort_values("pident_max", ascending=True)

    fig, ax = plt.subplots(figsize=(6, max(8, len(df) * 0.22)))
    colors = df["pident_max"].apply(
        lambda x: "#2ecc71" if x >= 90 else ("#f39c12" if x >= 75 else "#e74c3c"))
    ax.barh(df["strain"], df["pident_max"], color=colors, edgecolor="white",
            height=0.8)
    ax.axvline(90, color="green",  linestyle="--", alpha=0.5, linewidth=1,
               label="90% identity")
    ax.axvline(75, color="orange", linestyle="--", alpha=0.5, linewidth=1,
               label="75% identity")
    ax.set_xlabel("% nucleotide identity to query", fontsize=10)
    ax.set_ylabel("Strain", fontsize=10)
    ax.set_title(f"{gene} — BLASTn identity per strain\n(primary copy)",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(55, 102)
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    plt.tight_layout()
    out = FIGS_DIR / f"{gene}_identity_barplot.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating supplementary figures...\n")

    print("[Chromosomal locations]")
    for gene in GENES:
        plot_chromosomal_location(gene)

    print("\n[ML distance heatmaps]")
    plot_mldist_heatmap("KHS1", "KHS1_primary_iqtree")
    plot_mldist_heatmap("KHR1", "KHR1_iqtree")

    print("\n[Identity barplots]")
    for gene in GENES:
        plot_identity_barplot(gene)

    print(f"\nAll supplementary figures → {FIGS_DIR}")
