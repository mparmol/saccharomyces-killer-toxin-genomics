#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

PROJECT    = Path("")
SEQ_DIR    = PROJECT / "results" / "sequences"
META_TSV   = PROJECT / "results" / "genome_metadata.tsv"
SPP_CSV    = PROJECT / "Spp_Color_codes.csv"
DETECT_TSV = PROJECT / "results_2" / "discovery" / "detection_method_summary.tsv"
OUT_DIR    = PROJECT / "results_2" / "copy_number"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_meta = pd.read_csv(META_TSV, sep="\t")[["Strain", "Species"]].copy()
df_meta["Strain"] = df_meta["Strain"].str.strip()
df_meta["Species"] = df_meta["Species"].str.strip()
strain2spp = dict(zip(df_meta["Strain"], df_meta["Species"]))

df_spp = pd.read_csv(SPP_CSV)
spp2hex = dict(zip(df_spp["Species"].str.strip(), df_spp["HEX code"].str.strip()))
UNKNOWN_COLOR = "#95a5a6"

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
    "Saccharomyces westchinensis",
]

def analyze_gene(gene):
    detail_path = SEQ_DIR / gene / f"{gene}_blast_hits_detail.tsv"
    pa_path     = SEQ_DIR / gene / f"{gene}_presence_absence.tsv"
    if not detail_path.exists():
        print(f"  [SKIP] {gene}: {detail_path} no")
        return

    df = pd.read_csv(detail_path, sep="\t")
    df_pa = pd.read_csv(pa_path, sep="\t") if pa_path.exists() else None

    # Filter to valid strains only (removes KHR1 false positives)
    if df_pa is not None and "present" in df_pa.columns:
        valid_strains = set(df_pa.loc[df_pa["present"] == "present", "strain"])
        df = df[df["strain"].isin(valid_strains)]

    # Add nhmmer novel candidates from detection_method_summary.tsv
    if DETECT_TSV.exists():
        dt = pd.read_csv(DETECT_TSV, sep="\t")
        novel = dt[(dt["gene"] == gene) & (dt["is_novel"] == True)].copy()
        if len(novel) > 0:
            novel_rows = pd.DataFrame({
                "strain": novel["strain"].values,
                "chrom":  novel["chrom"].values,
                "copy":   1,
                "pident_max": np.nan,
                "detection_method": novel["detection_method"].values,
            })
            df = pd.concat([df, novel_rows], ignore_index=True)
            print(f"  + {len(novel)}")

    print(f"\n[{gene}] {len(df)} hits ({df['strain'].nunique()})")

    per_strain = (df.groupby("strain")
                    .agg(
                        n_copies    = ("copy", "max"),
                        chroms      = ("chrom", lambda x: ";".join(sorted(set(x)))),
                        n_chroms    = ("chrom", "nunique"),
                        pident_mean = ("pident_max", "mean"),
                        pident_min  = ("pident_max", "min"),
                        pident_max  = ("pident_max", "max"),
                    )
                    .reset_index())
    per_strain["species"] = per_strain["strain"].map(strain2spp)

    all_strains = df_meta["Strain"].tolist()
    present = set(per_strain["strain"])
    absent = [s for s in all_strains if s not in present]
    if absent:
        absent_rows = pd.DataFrame({
            "strain": absent,
            "n_copies": 0,
            "chroms": "",
            "n_chroms": 0,
            "pident_mean": np.nan,
            "pident_min": np.nan,
            "pident_max": np.nan,
            "species": [strain2spp.get(s, "Unknown") for s in absent],
        })
        per_strain = pd.concat([per_strain, absent_rows], ignore_index=True)

    per_strain["species"] = pd.Categorical(
        per_strain["species"], categories=SPP_ORDER, ordered=True
    )
    per_strain = per_strain.sort_values(["species", "strain"]).reset_index(drop=True)

    out_tsv = OUT_DIR / f"{gene}_copy_number_summary.tsv"
    per_strain.to_csv(out_tsv, sep="\t", index=False)
    print(f"  : {out_tsv}")

    df_hits = df.copy()
    pivot = (df_hits.groupby(["strain", "chrom"])
                    .size()
                    .reset_index(name="n_copies")
                    .pivot(index="strain", columns="chrom", values="n_copies")
                    .fillna(0)
                    .astype(int))

    def chr_sort_key(c):
        m = __import__("re").search(r"(\d+|I{1,3}|IV|V{1,3}|VI{1,3}|IX|X{1,3}|XI{1,3}|XIV|XV|XVI)$", c, __import__("re").I)
        roman = {"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8,
                 "IX":9,"X":10,"XI":11,"XII":12,"XIII":13,"XIV":14,"XV":15,"XVI":16}
        if m:
            v = m.group(1).upper()
            return roman.get(v, 99)
        return 99
    cols = sorted(pivot.columns, key=chr_sort_key)
    pivot = pivot[cols]

    pivot["species"] = pivot.index.map(strain2spp)
    pivot["species"] = pd.Categorical(pivot["species"], categories=SPP_ORDER, ordered=True)
    pivot = pivot.sort_values("species").drop(columns="species")

    n_strains = len(pivot)
    n_chroms  = len(pivot.columns)
    fig, ax = plt.subplots(figsize=(max(8, n_chroms * 0.7), max(12, n_strains * 0.28)))

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "copies", ["#f0f0f0", "#6baed6", "#08519c"], N=4
    )
    im = sns.heatmap(
        pivot, ax=ax, cmap=cmap, vmin=0, vmax=max(3, pivot.values.max()),
        linewidths=0.3, linecolor="white",
        annot=True, fmt="d", annot_kws={"size": 7},
        cbar_kws={"label": "Nº copias", "shrink": 0.4},
    )

    ylabels = ax.get_yticklabels()
    for lbl in ylabels:
        s = lbl.get_text()
        spp = strain2spp.get(s, "Unknown")
        lbl.set_color(spp2hex.get(spp, UNKNOWN_COLOR))
        lbl.set_fontsize(7)

    current_spp = None
    separator_pos = []
    strains_order = list(pivot.index)
    for i, s in enumerate(strains_order):
        spp = strain2spp.get(s, "Unknown")
        if spp != current_spp and current_spp is not None:
            separator_pos.append(i)
        current_spp = spp
    for pos in separator_pos:
        ax.axhline(pos, color="black", lw=1.2, alpha=0.6)

    ax.set_title(f"{gene} — ",
                 fontsize=12, fontweight="bold", pad=10)
    ax.set_xlabel("Cromosom", fontsize=10)
    ax.set_ylabel("Cepa", fontsize=10)
    ax.tick_params(axis="x", labelsize=8, rotation=45)

    spp_in_plot = per_strain[per_strain["n_copies"] > 0]["species"].dropna().unique()
    patches = [
        mpatches.Patch(color=spp2hex.get(s, UNKNOWN_COLOR),
                       label=s.replace("Saccharomyces ", "S. "))
        for s in SPP_ORDER if s in spp_in_plot
    ]
    ax.legend(handles=patches, bbox_to_anchor=(1.12, 1), loc="upper left",
              fontsize=7.5, title="Especie", title_fontsize=8,
              framealpha=0.9, edgecolor="gray")

    plt.tight_layout()
    out_pdf = OUT_DIR / f"{gene}_copy_heatmap.pdf"
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.savefig(str(out_pdf).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Heatmap: {out_pdf}")

    copy_counts = per_strain["n_copies"].value_counts().sort_index()
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    bars = ax2.bar(copy_counts.index.astype(str), copy_counts.values,
                   color=["#f0f0f0", "#6baed6", "#08519c", "#08306b"][:len(copy_counts)],
                   edgecolor="white", linewidth=0.8)
    for bar, v in zip(bars, copy_counts.values):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                 str(v), ha="center", va="bottom", fontsize=9)
    ax2.set_xlabel("Nm", fontsize=10)
    ax2.set_ylabel("Nm", fontsize=10)
    ax2.set_title(f"{gene} — ",
                  fontsize=11, fontweight="bold")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    out_bar = OUT_DIR / f"{gene}_copy_barplot.pdf"
    plt.savefig(out_bar, dpi=300, bbox_inches="tight")
    plt.savefig(str(out_bar).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Barplot: {out_bar}")

    return per_strain


def combined_table(results):
    rows = []
    for gene, df_s in results.items():
        for _, r in df_s.iterrows():
            rows.append({
                "Strain":   r["strain"],
                "Species":  r["species"],
                f"{gene}_copies": int(r["n_copies"]),
                f"{gene}_chroms": r["chroms"],
            })
    df_khs1 = results.get("KHS1")
    df_khr1 = results.get("KHR1")
    if df_khs1 is not None and df_khr1 is not None:
        df_khs1 = df_khs1[["strain", "species", "n_copies", "chroms"]].rename(
            columns={"n_copies": "KHS1_copies", "chroms": "KHS1_chroms"})
        df_khr1 = df_khr1[["strain", "n_copies", "chroms"]].rename(
            columns={"n_copies": "KHR1_copies", "chroms": "KHR1_chroms"})
        df_comb = df_khs1.merge(df_khr1, on="strain", how="outer")
        df_comb["species"] = pd.Categorical(df_comb["species"], categories=SPP_ORDER, ordered=True)
        df_comb = df_comb.sort_values(["species", "strain"])
        out = OUT_DIR / "copy_number_combined.tsv"
        df_comb.to_csv(out, sep="\t", index=False)
        print(f"\nTabl: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Copy number analysis — KHS1 y KHR1")
    print("=" * 60)

    results = {}
    for gene in ["KHS1", "KHR1"]:
        res = analyze_gene(gene)
        if res is not None:
            results[gene] = res

    combined_table(results)
    print(f"\nOutput: {OUT_DIR}")
