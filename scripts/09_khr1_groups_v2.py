#!/usr/bin/env python3
"""
KHR1 group analysis v2 — 68 clean sequences (C1 + C2 + minor singletons).

Uses the full 68-seq tree from KHR1_iqtree.contree (built by 10_fix_khr1_full.sh).
Groups: C1 = chrIX (canonical, SK1-like), C2 = chrIII, singletons = chrXV/chrXIV.
"""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

PROJECT  = Path("")
SEQ_DIR  = PROJECT / "results" / "sequences"  / "KHR1"
TREE_DIR = PROJECT / "results" / "trees"      / "KHR1"
FIGS_DIR = PROJECT / "results" / "figures"
FIGS_DIR.mkdir(exist_ok=True)

VALID_CHROMS = {"chrIX", "chrIII", "chrXV", "chrXIV"}
GROUP_MAP    = {"chrIX": "C1", "chrIII": "C2", "chrXV": "C3", "chrXIV": "C3"}
COLOR_MAP    = {"C1": "#2980b9", "C2": "#e74c3c", "C3": "#27ae60", "other": "#95a5a6"}

# ── Load strain → chrom / group mapping ──────────────────────────────────────
df = pd.read_csv(SEQ_DIR / "KHR1_blast_hits_detail.tsv", sep="\t")
df = df[df["copy"] == 1].copy()
df = df[df["chrom"].isin(VALID_CHROMS)]

chrom_of  = dict(zip(df["strain"], df["chrom"]))
group_of  = {s: GROUP_MAP[c] for s, c in chrom_of.items()}

print("=== KHR1 group composition (68 clean seqs) ===")
for grp in ["C1", "C2", "C3"]:
    members = [s for s, g in group_of.items() if g == grp]
    chroms  = sorted({chrom_of[s] for s in members})
    print(f"  {grp} ({', '.join(chroms)}): {len(members)} strains")
print()

# ── Read tree ────────────────────────────────────────────────────────────────
treefile = TREE_DIR / "KHR1_iqtree.contree"
if not treefile.exists():
    raise FileNotFoundError(f"Tree not found: {treefile}\nRun 10_fix_khr1_full.sh first.")

iqtree_rpt = TREE_DIR / "KHR1_iqtree.iqtree"
model_label = "HKY+F+R2"
if iqtree_rpt.exists():
    import re
    text = iqtree_rpt.read_text()
    m = re.search(r"Best-fit model according to BIC:\s+(\S+)", text)
    if m:
        model_label = m.group(1)

tree = Phylo.read(treefile, "newick")
tree.ladderize()

# ── Verify monophyly of C1 and C2 ────────────────────────────────────────────
tree2 = Phylo.read(treefile, "newick")
tree2.ladderize()

c1_tips = [t.name for t in tree2.get_terminals() if group_of.get(t.name, "") == "C1"]
c2_tips = [t.name for t in tree2.get_terminals() if group_of.get(t.name, "") == "C2"]

def check_monophyly(tree, tips, group_name):
    if len(tips) < 2:
        print(f"  {group_name}: only {len(tips)} tip(s), skipping monophyly check")
        return
    mrca      = tree.common_ancestor(tips)
    sub_tips  = [t.name for t in mrca.get_terminals()]
    intruders = [s for s in sub_tips if group_of.get(s, "") != group_name]
    pure      = len(intruders) == 0
    conf      = mrca.confidence
    try:
        bs = int(float(str(conf).split("/")[0])) if conf else None
    except Exception:
        bs = None
    print(f"  {group_name}: {len(tips)} tips, monophyletic={pure}, root_UFBoot={bs}")
    if intruders:
        print(f"    Intruders: {intruders}")

print("=== Monophyly check ===")
check_monophyly(tree2, c1_tips, "C1")
check_monophyly(tree2, c2_tips, "C2")
print()

# ── Draw annotated tree ──────────────────────────────────────────────────────
def tip_color(name):
    if not name:
        return "black"
    s = name.lstrip("_R_")
    return COLOR_MAP.get(group_of.get(s, "other"), "#95a5a6")

n_tips = tree.count_terminals()
fig, ax = plt.subplots(figsize=(14, max(14, n_tips * 0.25)))

Phylo.draw(tree, axes=ax, do_show=False,
           label_func=lambda c: c.name if c.is_terminal() else "",
           label_colors=lambda name: tip_color(name))

ax.set_xlabel("Substitutions per site", fontsize=10)
ax.set_title(
    f"KHR1 killer toxin — {n_tips} Saccharomyces strains\n"
    "C1: chrIX  ·  C2: chrIII  ·  C3: singletons (chrXV, chrXIV)",
    fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.yaxis.set_visible(False)

c1_n = len([s for s, g in group_of.items() if g == "C1"])
c2_n = len([s for s, g in group_of.items() if g == "C2"])
c3_n = len([s for s, g in group_of.items() if g == "C3"])

patches = [
    mpatches.Patch(color=COLOR_MAP["C1"], label=f"C1 – chrIX, n={c1_n} (SK1-like)"),
    mpatches.Patch(color=COLOR_MAP["C2"], label=f"C2 – chrIII, n={c2_n}"),
    mpatches.Patch(color=COLOR_MAP["C3"], label=f"C3 – chrXV/XIV, n={c3_n} (singletons)"),
]
ax.legend(handles=patches, loc="lower left", fontsize=10,
          title="KHR1 group", title_fontsize=9)

ax.text(0.98, 0.98,
        f"Model: {model_label}\nUFBoot 1000 | SH-aLRT 1000 | IQ-TREE 3.1.2",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.9))

plt.tight_layout()
out_pdf = FIGS_DIR / "KHR1_C1C2_tree_annotated.pdf"
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.savefig(str(out_pdf).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Annotated tree → {out_pdf}")

# ── Summary table ────────────────────────────────────────────────────────────
rows = []
for grp, chrom_label, notes in [
    ("C1", "chrIX",       "Canonical (SK1-like)"),
    ("C2", "chrIII",      "Likely translocation"),
    ("C3", "chrXV, chrXIV", "Singletons, marginal"),
]:
    members = df[df["chrom"].isin(
        {c for c, g in GROUP_MAP.items() if g == grp}
    )]
    if members.empty:
        continue
    rows.append({
        "Group":         grp,
        "Chromosome":    chrom_label,
        "N strains":     len(members),
        "% identity (mean)": f"{members['pident_max'].mean():.1f}",
        "% identity (range)": f"{members['pident_max'].min():.1f}–{members['pident_max'].max():.1f}",
        "Notes":         notes,
    })

# Excluded (false positives — recorded for transparency)
excluded = pd.read_csv(SEQ_DIR / "KHR1_blast_hits_detail.tsv", sep="\t")
excluded = excluded[(excluded["copy"] == 1) & (~excluded["chrom"].isin(VALID_CHROMS))]
if not excluded.empty:
    rows.append({
        "Group":         "Excluded",
        "Chromosome":    ", ".join(sorted(excluded["chrom"].unique())),
        "N strains":     len(excluded),
        "% identity (mean)": f"{excluded['pident_max'].mean():.1f}",
        "% identity (range)": f"{excluded['pident_max'].min():.1f}–{excluded['pident_max'].max():.1f}",
        "Notes":         "tBLASTx-only, <41% aa identity — false positives",
    })

df_g = pd.DataFrame(rows)
print("\n=== KHR1 group summary ===")
print(df_g.to_string(index=False))

df_g.to_csv(FIGS_DIR / "KHR1_group_summary.tsv",  sep="\t", index=False)
df_g.to_excel(FIGS_DIR / "KHR1_group_summary.xlsx", index=False)

# Figure table
fig, ax = plt.subplots(figsize=(13, 1.5 + 0.6 * len(rows)))
ax.axis("off")
tbl = ax.table(
    cellText=df_g.values.tolist(),
    colLabels=list(df_g.columns),
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
ax.set_title("KHR1 group summary", fontsize=11, fontweight="bold", pad=10)
plt.tight_layout()
out_tbl = FIGS_DIR / "KHR1_group_summary_table.pdf"
plt.savefig(out_tbl, dpi=300, bbox_inches="tight")
plt.savefig(str(out_tbl).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Summary table → {out_tbl}")
