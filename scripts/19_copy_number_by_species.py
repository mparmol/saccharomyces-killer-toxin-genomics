
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

PROJ    = Path(r"")
META    = PROJ / "" / "01_metadata" / "genome_metadata.tsv"
CN_DIR  = PROJ / "" / ""
SPP_CSV = PROJ / "Spp_Color_codes.csv"
OUT     = PROJ / "" / "28-5-26_2" / "" / "" / ""
OUT.mkdir(parents=True, exist_ok=True)

# ── species colours (same palette as other figures) ────────────────────────
_sdf = pd.read_csv(SPP_CSV)
SPP2HEX = dict(zip(_sdf["Species"].str.strip(), _sdf["HEX code"].str.strip()))

# ── species order ──────────────────────────────────────────────────────────
# Phylogenetic order requested by reviewer (DPN, comment #138):
# Scer > Spar > Smik > Sjur > Skud > Sarb > Suva > Schi > Seub
SPECIES_ORDER = [
    "Saccharomyces cerevisiae",
    "Saccharomyces paradoxus",
    "Saccharomyces mikatae",
    "Saccharomyces jurei",
    "Saccharomyces kudriavzevii",
    "Saccharomyces arboricola",
    "Saccharomyces uvarum",
    "Saccharomyces chiloensis",
    "Saccharomyces eubayanus",
]
SPP_LABEL = {s: f"S. {s.split()[1]}" for s in SPECIES_ORDER}
N_STRAINS = {                    # total strains per species
    "Saccharomyces arboricola":   8,
    "Saccharomyces cerevisiae":  21,
    "Saccharomyces chiloensis":   2,
    "Saccharomyces eubayanus":   11,
    "Saccharomyces jurei":        2,
    "Saccharomyces kudriavzevii": 7,
    "Saccharomyces mikatae":     10,
    "Saccharomyces paradoxus":   15,
    "Saccharomyces uvarum":       7,
}

# ── chromosomes ────────────────────────────────────────────────────────────
ROMAN = ["I","II","III","IV","V","VI","VII","VIII",
         "IX","X","XI","XII","XIII","XIV","XV","XVI"]
CHR_ORDER = [f"chr{r}" for r in ROMAN]

# ── strain -> species ──────────────────────────────────────────────────────
meta = pd.read_csv(META, sep="\t")[["Strain","Species"]].copy()
meta["Strain"] = meta["Strain"].str.strip()
CORR = {
    "UFRJ50816T":"UFRJ50816","CDFM21L.1":"CDFM21L1","MNFM4L.2":"MNFM4L2",
    "RTSP5L.1":"RTSP5L1","SN-QL4-2":"SNQL42","LL2011-012":"LL2011012",
    "LL2012-001":"LL2012001","LL2012-016":"LL2012016","LL2012-018":"LL2012018",
    "MSH-604":"MSH604","BR6-2":"BR62","JXXY16.1":"JXXY161",
    "UWOPS03-461.4":"UWOPS034614","UWOPS91-917.1":"UWOPS919171","XXYS1.4":"XXYS14",
    "EM14S01-3B":"EM14S013B",
}
for old, new in CORR.items():
    meta["Strain"] = meta["Strain"].replace(old, new)
s2spp = dict(zip(meta["Strain"], meta["Species"]))


def build_matrices(gene):
    """
    Returns two DataFrames (index=CHR_ORDER, columns=SPECIES_ORDER):
      mat_copies : total copies per species x chromosome
      mat_strains: number of strains with >= 1 copy on that chromosome
    """
    df = pd.read_csv(CN_DIR / f"{gene}_copy_number_summary.tsv", sep="\t")
    df["strain"] = df["strain"].str.strip()
    if "species" in df.columns:
        df["Species"] = df["species"].str.strip()
        mask = df["Species"].isna() | (df["Species"] == "")
        df.loc[mask, "Species"] = df.loc[mask, "strain"].map(s2spp)
    else:
        df["Species"] = df["strain"].map(s2spp)

    mat_copies  = pd.DataFrame(0, index=CHR_ORDER, columns=SPECIES_ORDER, dtype=float)
    mat_strains = pd.DataFrame(0, index=CHR_ORDER, columns=SPECIES_ORDER, dtype=int)

    for _, row in df.iterrows():
        spp = row["Species"]
        if spp not in SPECIES_ORDER:
            continue
        chroms_raw = str(row.get("chroms","")).strip()
        if not chroms_raw or chroms_raw == "nan":
            continue
        n_copies  = int(row["n_copies"])
        n_chroms  = max(int(row["n_chroms"]), 1)
        chrom_list = [c.strip() for c in chroms_raw.split(";") if c.strip() in CHR_ORDER]
        for ch in chrom_list:
            # copies on this specific chromosome:
            # if multiple chromosomes -> 1 copy per chromosome (translocation)
            # if single chromosome   -> all copies are there (dimeric locus)
            copies_here = n_copies if n_chroms == 1 else 1
            mat_copies.loc[ch, spp]  += copies_here
            mat_strains.loc[ch, spp] += 1

    return mat_copies, mat_strains


# ── colormaps ──────────────────────────────────────────────────────────────
cmap_khs1 = LinearSegmentedColormap.from_list(
    "khs1", ["#ffffff","#fddbc7","#f4a582","#d6604d","#b2182b"], N=256)
cmap_khr1 = LinearSegmentedColormap.from_list(
    "khr1", ["#ffffff","#d1e5f0","#92c5de","#4393c3","#2166ac"], N=256)
CMAPS  = {"KHS1": cmap_khs1, "KHR1": cmap_khr1}
TITLES = {"KHS1": "khs1",    "KHR1": "khr1"}

# ── figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 7),
                         gridspec_kw={"wspace": 0.55})
fig.patch.set_facecolor("white")

PANEL = {"KHS1": "A", "KHR1": "B"}
for ax, gene in zip(axes, ("KHS1","KHR1")):
    mat_c, mat_s = build_matrices(gene)
    vmax = mat_c.values.max() if mat_c.values.max() > 0 else 1

    im = ax.imshow(mat_c.values, aspect="auto", cmap=CMAPS[gene],
                   vmin=0, vmax=vmax, interpolation="nearest")

    # minor grid lines between cells
    ax.set_xticks(np.arange(-0.5, len(SPECIES_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CHR_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    # axes — species labels coloured by species, italic
    ax.set_xticks(range(len(SPECIES_ORDER)))
    ax.set_xticklabels([SPP_LABEL[s] for s in SPECIES_ORDER],
                       rotation=45, ha="right", fontstyle="italic", fontsize=9.5)
    for tick, spp in zip(ax.get_xticklabels(), SPECIES_ORDER):
        tick.set_color(SPP2HEX.get(spp, "#555555"))
    ax.set_yticks(range(len(CHR_ORDER)))
    ax.set_yticklabels([f"chr{r}" for r in ROMAN], fontsize=9)

    # cell annotations: "N strains / total  (X copies)" on two lines
    for r, ch in enumerate(CHR_ORDER):
        for c, spp in enumerate(SPECIES_ORDER):
            n_copies = mat_c.iloc[r, c]
            n_strains = mat_s.iloc[r, c]
            if n_copies == 0:
                continue
            total = N_STRAINS[spp]
            brightness = n_copies / vmax
            fc = "white" if brightness > 0.55 else "#222222"

            # top line: strains
            strain_lbl = f"{n_strains}/{total}"
            # bottom line: copies (only show if different from strains)
            if int(n_copies) != n_strains:
                copy_lbl = f"({int(n_copies)} cop.)"
                lbl = f"{strain_lbl}\n{copy_lbl}"
                fs_top, fs_bot = 8, 6.5
            else:
                lbl = strain_lbl
                fs_top = 8.5
                fs_bot = None

            if fs_bot:
                # two-line text
                ax.text(c, r - 0.12, strain_lbl,
                        ha="center", va="center", fontsize=fs_top,
                        color=fc, fontweight="bold")
                ax.text(c, r + 0.28, copy_lbl,
                        ha="center", va="center", fontsize=fs_bot,
                        color=fc)
            else:
                ax.text(c, r, lbl,
                        ha="center", va="center", fontsize=fs_top,
                        color=fc, fontweight="bold")

    # colorbar
    cb = fig.colorbar(im, ax=ax, shrink=0.52, pad=0.02)
    cb.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    cb.ax.tick_params(labelsize=8)
    cb.set_label("Total copies", fontsize=9)

    # title -- uppercase italic gene name
    ax.set_title(f"$\\mathit{{{gene}}}$", fontsize=14,
                 fontweight="bold", pad=10)
    # panel letter (A / B) requested by reviewer (#170)
    ax.text(-0.18, 1.04, PANEL[gene], transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="bottom", ha="left")
    ax.spines[:].set_visible(False)

# shared legend note
fig.text(0.5, -0.03,
         "Cell text: N strains / species total  (copies shown in parentheses when > strains)",
         ha="center", fontsize=8.5, color="#555555")

for fmt in ("pdf","png","svg"):
    fig.savefig(OUT / f"copy_number_species_chrom.{fmt}",
                dpi=220, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("Done.")
