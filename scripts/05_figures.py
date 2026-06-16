#!/usr/bin/env python3
"""
Step 5: Publication-quality figures and summary tables.

Generates:
  1. Presence/absence heatmap for KHS1 + KHR1 across all 83 strains
  2. BLAST identity distribution plots per gene
  3. Phylogenetic tree visualizations (requires ete3 or uses ASCII output)
  4. Alignment statistics table (length, variable sites, informative sites)
  5. Combined summary table (strain info + presence/absence + identity)
"""

from pathlib import Path
import csv
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import seaborn as sns
from Bio import SeqIO, AlignIO

PROJECT_DIR  = Path("")
RESULTS_DIR  = PROJECT_DIR / "results"
SEQ_DIR      = RESULTS_DIR / "sequences"
ALN_DIR      = RESULTS_DIR / "alignments"
TREE_DIR     = RESULTS_DIR / "trees"
FIGS_DIR     = RESULTS_DIR / "figures"
FIGS_DIR.mkdir(exist_ok=True)

GENES = ["KHS1", "KHR1"]

# ── colour palette ───────────────────────────────────────────────────────────
PALETTE = {
    "present": "#2ecc71",
    "absent":  "#e74c3c",
    "bg":      "#f8f9fa",
}


# ════════════════════════════════════════════════════════════════════════════
# 1. Load presence/absence data
# ════════════════════════════════════════════════════════════════════════════
def load_pa(gene):
    pa_file = SEQ_DIR / gene / f"{gene}_presence_absence.tsv"
    if not pa_file.exists():
        return pd.DataFrame()
    return pd.read_csv(pa_file, sep="\t")


def build_pa_matrix():
    dfs = {}
    for gene in GENES:
        df = load_pa(gene)
        if df.empty:
            continue
        dfs[gene] = df.set_index("strain")["present"]

    if not dfs:
        return pd.DataFrame()

    mat = pd.DataFrame(dfs)
    mat = mat.applymap(lambda x: 1 if x == "present" else 0)
    return mat


# ════════════════════════════════════════════════════════════════════════════
# 2. Presence/absence heatmap
# ════════════════════════════════════════════════════════════════════════════
def plot_presence_absence(mat):
    if mat.empty:
        return

    fig, ax = plt.subplots(figsize=(4, max(8, len(mat) * 0.22)))

    cmap = ListedColormap([PALETTE["absent"], PALETTE["present"]])
    sns.heatmap(
        mat, ax=ax, cmap=cmap, linewidths=0.5, linecolor="white",
        cbar=False, vmin=0, vmax=1,
        xticklabels=True, yticklabels=True,
    )

    ax.set_title("Killer toxin presence/absence\nacross 83 Saccharomyces strains",
                 fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=10, rotation=0)
    ax.tick_params(axis="y", labelsize=7.5)

    present_patch = mpatches.Patch(color=PALETTE["present"], label="Present")
    absent_patch  = mpatches.Patch(color=PALETTE["absent"],  label="Absent")
    ax.legend(handles=[present_patch, absent_patch],
              loc="lower right", bbox_to_anchor=(1.22, 0), fontsize=9)

    plt.tight_layout()
    out = FIGS_DIR / "presence_absence_heatmap.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════════
# 3. BLAST identity violin/box plots
# ════════════════════════════════════════════════════════════════════════════
def plot_identity_distribution():
    fig, axes = plt.subplots(1, len(GENES), figsize=(4 * len(GENES), 5))
    if len(GENES) == 1:
        axes = [axes]

    for ax, gene in zip(axes, GENES):
        detail_file = SEQ_DIR / gene / f"{gene}_blast_hits_detail.tsv"
        pa_file     = SEQ_DIR / gene / f"{gene}_presence_absence.tsv"
        if not detail_file.exists():
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        df = pd.read_csv(detail_file, sep="\t")
        # Filter to strains confirmed present (removes tBLASTx false positives)
        if pa_file.exists():
            pa = pd.read_csv(pa_file, sep="\t")
            valid = set(pa[pa["present"] == "present"]["strain"])
            df = df[df["strain"].isin(valid)]
        ax.violinplot(df["pident_max"], positions=[0], showmedians=True,
                      widths=0.6)
        ax.boxplot(df["pident_max"], positions=[0], widths=0.15,
                   patch_artist=True,
                   boxprops=dict(facecolor="white", linewidth=1.5),
                   medianprops=dict(color="red", linewidth=2),
                   whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=1.5),
                   flierprops=dict(marker="o", markersize=3, alpha=0.5))

        ax.set_title(f"{gene}", fontsize=11, fontweight="bold")
        ax.set_ylabel("% identity to query", fontsize=9)
        ax.set_ylim(50, 102)
        ax.set_xticks([])
        ax.yaxis.grid(True, linestyle="--", alpha=0.5)
        ax.set_facecolor(PALETTE["bg"])

        # Annotate stats
        med = df["pident_max"].median()
        mn  = df["pident_max"].min()
        mx  = df["pident_max"].max()
        ax.text(0.97, 0.97, f"n={len(df)}\nmed={med:.1f}%\nmin={mn:.1f}%\nmax={mx:.1f}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.suptitle("BLASTn % identity distribution", fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = FIGS_DIR / "blast_identity_distribution.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════════
# 4. Alignment statistics
# ════════════════════════════════════════════════════════════════════════════
ALN_TRIM_NAME = {
    "KHS1": "KHS1_mafft_trim.fasta",
    "KHR1": "KHR1_mafft_gt50.fasta",
}

def alignment_stats(gene):
    aln_file = ALN_DIR / gene / ALN_TRIM_NAME.get(gene, f"{gene}_mafft_trim.fasta")
    if not aln_file.exists():
        return None

    aln = AlignIO.read(aln_file, "fasta")
    n_seq   = len(aln)
    aln_len = aln.get_alignment_length()

    # Variable sites
    variable = 0
    parsimony = 0
    for i in range(aln_len):
        col = aln[:, i]
        bases = set(col.upper().replace("-", ""))
        if len(bases) > 1:
            variable += 1
            counts = {b: col.upper().count(b) for b in bases}
            if sum(1 for c in counts.values() if c >= 2) >= 2:
                parsimony += 1

    return {
        "Gene": gene,
        "N sequences": n_seq,
        "Alignment length (bp)": aln_len,
        "Variable sites": variable,
        "Parsimony-informative sites": parsimony,
        "% variable": f"{100*variable/aln_len:.1f}",
        "% parsimony-informative": f"{100*parsimony/aln_len:.1f}",
    }


def write_alignment_stats_table():
    rows = [r for r in (alignment_stats(g) for g in GENES) if r is not None]
    if not rows:
        return

    df = pd.DataFrame(rows)
    out_tsv = FIGS_DIR / "alignment_statistics.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Saved: {out_tsv}")

    # Also render as a figure/table
    fig, ax = plt.subplots(figsize=(9, 1.5 + 0.4 * len(rows)))
    ax.axis("off")
    tbl = ax.table(
        cellText=df.values, colLabels=df.columns,
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#ecf0f1")
    plt.title("Alignment Statistics", fontsize=11, fontweight="bold", pad=8)
    plt.tight_layout()
    out_fig = FIGS_DIR / "alignment_statistics_table.pdf"
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.savefig(str(out_fig).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_fig}")


# ════════════════════════════════════════════════════════════════════════════
# 5. Combined summary table
# ════════════════════════════════════════════════════════════════════════════
def write_combined_summary():
    mat = build_pa_matrix()
    if mat.empty:
        return

    rows = []
    for strain in mat.index:
        row = {"Strain": strain}
        for gene in GENES:
            pa_df = load_pa(gene)
            if pa_df.empty:
                row[f"{gene}_presence"] = "NA"
                row[f"{gene}_pident"]   = "NA"
                row[f"{gene}_copies"]   = "NA"
                continue
            pa_df = pa_df.set_index("strain")
            if strain in pa_df.index:
                r = pa_df.loc[strain]
                row[f"{gene}_presence"] = r["present"]
                row[f"{gene}_pident"]   = r["best_pident"]
                row[f"{gene}_copies"]   = r["n_copies"]
            else:
                row[f"{gene}_presence"] = "absent"
                row[f"{gene}_pident"]   = "NA"
                row[f"{gene}_copies"]   = 0
        rows.append(row)

    df = pd.DataFrame(rows)
    out = FIGS_DIR / "combined_summary_table.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"  Saved: {out}")

    # Excel-friendly version
    df.to_excel(str(out).replace(".tsv", ".xlsx"), index=False)
    print(f"  Saved: {str(out).replace('.tsv', '.xlsx')}")


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating figures and tables...")

    mat = build_pa_matrix()
    print("\n[1] Presence/absence heatmap")
    plot_presence_absence(mat)

    print("\n[2] BLAST identity distributions")
    plot_identity_distribution()

    print("\n[3] Alignment statistics")
    write_alignment_stats_table()

    print("\n[4] Combined summary table")
    write_combined_summary()

    print(f"\nAll figures → {FIGS_DIR}")
