
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

PROJ    = Path(r"")
DET_DIR = PROJ / "results" / "02_blast_detection"
META    = PROJ / "results" / "01_metadata" / "genome_metadata.tsv"
SPP_CSV = PROJ / "Spp_Color_codes.csv"
OUT     = PROJ / "" / ""
OUT.mkdir(exist_ok=True)

# ── species colour map ─────────────────────────────────────────────────────
spp_df  = pd.read_csv(SPP_CSV)
spp2hex = dict(zip(spp_df["Species"].str.strip(),
                   spp_df["HEX code"].str.strip()))
UNKNOWN = "#aaaaaa"

SPP_ORDER = [
    "Saccharomyces cerevisiae",
    "Saccharomyces paradoxus",
    "Saccharomyces mikatae",
    "Saccharomyces jurei",
    "Saccharomyces kudriavzevii",
    "Saccharomyces arboricola",
    "Saccharomyces eubayanus",
    "Saccharomyces uvarum",
    "Saccharomyces chiloensis",
]

# ── strain -> species (with name corrections) ──────────────────────────────
CORR = {
    "UFRJ50816T": "UFRJ50816",
    "CDFM21L.1":  "CDFM21L1",
    "MNFM4L.2":   "MNFM4L2",
    "RTSP5L.1":   "RTSP5L1",
    "SN-QL4-2":   "SNQL42",
    "LL2011-012":  "LL2011012",
    "LL2012-001":  "LL2012001",
    "LL2012-016":  "LL2012016",
    "LL2012-018":  "LL2012018",
    "MSH-604":    "MSH604",
    "BR6-2":      "BR62",
    "JXXY16.1":   "JXXY161",
    "UWOPS03-461.4": "UWOPS034614",
    "UWOPS91-917.1": "UWOPS919171",
    "XXYS1.4":    "XXYS14",
    "EM14S01-3B": "EM14S013B",
}
meta = pd.read_csv(META, sep="\t")[["Strain","Species"]].copy()
meta["Strain"] = meta["Strain"].str.strip().replace(CORR)
meta["Species"] = meta["Species"].str.strip()
s2spp = dict(zip(meta["Strain"], meta["Species"]))      # normalised -> species
all_strains_norm = meta["Strain"].tolist()

# ── chromosome order ───────────────────────────────────────────────────────
ROMAN = ["I","II","III","IV","V","VI","VII","VIII",
         "IX","X","XI","XII","XIII","XIV","XV","XVI"]
ROM_RANK = {f"chr{r}": i for i, r in enumerate(ROMAN)}

def chr_rank(c):
    return ROM_RANK.get(c, 99)


def make_heatmap(gene):
    detail = DET_DIR / gene / f"{gene}_blast_hits_detail.tsv"
    pa     = DET_DIR / gene / f"{gene}_presence_absence.tsv"

    df = pd.read_csv(detail, sep="\t")
    df["strain"] = df["strain"].str.strip()

    # keep only validated strains
    if pa.exists():
        df_pa = pd.read_csv(pa, sep="\t")
        valid = set(df_pa.loc[df_pa["present"] == "present", "strain"].str.strip())
        df = df[df["strain"].isin(valid)]

    # pivot: strains x chromosomes, value = number of hits (copies)
    pivot = (df.groupby(["strain","chrom"])
               .size()
               .reset_index(name="n")
               .pivot(index="strain", columns="chrom", values="n")
               .fillna(0).astype(int))

    # sort chromosomes
    ordered_cols = sorted(pivot.columns, key=chr_rank)
    pivot = pivot[ordered_cols]

    # add absent strains as all-zero rows
    present = set(pivot.index)
    absent  = [s for s in all_strains_norm if s not in present]
    if absent:
        zero = pd.DataFrame(0, index=absent, columns=pivot.columns)
        pivot = pd.concat([pivot, zero])

    # assign & sort by species
    pivot["_spp"] = pd.Categorical(
        pivot.index.map(s2spp),
        categories=SPP_ORDER, ordered=True
    )
    pivot = pivot.sort_values("_spp").drop(columns="_spp")

    # ── figure ──────────────────────────────────────────────────────────────
    n_s = len(pivot)
    n_c = len(pivot.columns)
    fig, ax = plt.subplots(figsize=(max(8, n_c * 0.75), max(12, n_s * 0.29)))
    fig.patch.set_facecolor("white")

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "copies", ["#f7fbff","#6baed6","#08519c"], N=256)
    vmax = max(2, pivot.values.max())

    im = ax.imshow(pivot.values, aspect="auto", cmap=cmap,
                   vmin=0, vmax=vmax, interpolation="nearest")

    # grid
    ax.set_xticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_s, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # cell annotations (only non-zero)
    for r in range(n_s):
        for c in range(n_c):
            v = pivot.values[r, c]
            if v > 0:
                fc = "white" if v / vmax > 0.55 else "#222222"
                ax.text(c, r, str(v), ha="center", va="center",
                        fontsize=6.5, color=fc, fontweight="bold")

    # x-axis (chromosomes)
    ax.set_xticks(range(n_c))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)

    # y-axis (strain labels coloured by species -- THE FIX)
    ax.set_yticks(range(n_s))
    ax.set_yticklabels(pivot.index, fontsize=6.5)
    for lbl in ax.get_yticklabels():
        strain = lbl.get_text()
        spp    = s2spp.get(strain, "Unknown")
        lbl.set_color(spp2hex.get(spp, UNKNOWN))

    # species separator lines
    prev_spp = None
    for i, strain in enumerate(pivot.index):
        spp = s2spp.get(strain, "Unknown")
        if prev_spp is not None and spp != prev_spp:
            ax.axhline(i - 0.5, color="black", lw=1.1, alpha=0.55)
        prev_spp = spp

    # colourbar
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.01)
    cb.set_label("No. copies", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # species legend
    spp_present = [s for s in SPP_ORDER
                   if any(s2spp.get(st) == s for st in pivot.index
                          if pivot.loc[st].sum() > 0)]
    patches = [mpatches.Patch(color=spp2hex.get(s, UNKNOWN),
                              label=f"S. {s.split()[1]}")
               for s in spp_present]
    ax.legend(handles=patches, bbox_to_anchor=(1.13, 1), loc="upper left",
              title="Species", title_fontsize=8.5,
              framealpha=0.9, edgecolor="#cccccc",
              prop={"size": 7.5, "style": "italic"})

    ax.set_title(f"$\\it{{{gene.lower()}}}$ — copy number per strain and chromosome",
                 fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel("Chromosome", fontsize=9)
    ax.set_ylabel("Strain", fontsize=9)
    ax.spines[:].set_visible(False)

    fig.tight_layout()
    for fmt in ("pdf","png","svg"):
        fig.savefig(OUT / f"{gene}_copy_heatmap_corrected.{fmt}",
                    dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"{gene}: saved.")


for gene in ("KHS1","KHR1"):
    make_heatmap(gene)
print("Done.")
