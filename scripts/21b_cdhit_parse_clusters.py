
import re
from pathlib import Path
import pandas as pd

PROJ = Path(r"")
CDHIT_DIR = PROJ / "res" / "" / "cdhit"
BLAST_DIR = PROJ / "res" / ""
META      = PROJ / "res" / "01_metadata" / "genome_metadata.tsv"

CORR = {
    "UFRJ50816T":"UFRJ50816","CDFM21L.1":"CDFM21L1","MNFM4L.2":"MNFM4L2",
    "RTSP5L.1":"RTSP5L1","SN-QL4-2":"SNQL42","LL2011-012":"LL2011012",
    "LL2012-001":"LL2012001","LL2012-016":"LL2012016","LL2012-018":"LL2012018",
    "MSH-604":"MSH604","BR6-2":"BR62","JXXY16.1":"JXXY161",
    "UWOPS03-461.4":"UWOPS034614","UWOPS91-917.1":"UWOPS919171","XXYS1.4":"XXYS14",
    "EM14S01-3B":"EM14S013B",
}

meta = pd.read_csv(META, sep="\t")[["Strain","Species"]].copy()
meta["Strain"] = meta["Strain"].str.strip().replace(CORR)
s2spp = dict(zip(meta["Strain"], meta["Species"]))


def load_strain_chrom(gene):
    """Return dict {strain: primary_chrom} from blast_hits_detail.tsv."""
    detail = BLAST_DIR / gene / f"{gene}_blast_hits_detail.tsv"
    pa     = BLAST_DIR / gene / f"{gene}_presence_absence.tsv"
    df = pd.read_csv(detail, sep="\t")
    df["strain"] = df["strain"].str.strip().replace(CORR)
    # Filter to valid (present) strains
    if pa.exists():
        dfp = pd.read_csv(pa, sep="\t")
        valid = set(dfp.loc[dfp["present"]=="present","strain"].str.strip().replace(CORR))
        df = df[df["strain"].isin(valid)]
    # Primary chrom = highest bitscore per strain
    df_best = df.sort_values("bitscore_max", ascending=False).drop_duplicates("strain")
    return dict(zip(df_best["strain"], df_best["chrom"]))


def parse_clstr(clstr_path):
    clusters = []
    current = None
    for line in open(clstr_path):
        line = line.rstrip()
        if line.startswith(">Cluster"):
            if current is not None:
                clusters.append(current)
            current = {"id": int(line.split()[1]), "representative": None, "members": []}
        elif current is not None and line:
            m = re.search(r'>([^.]+)\.\.\.', line)
            if m:
                seq_id = m.group(1)
                current["members"].append(seq_id)
                if "*" in line:
                    current["representative"] = seq_id
    if current is not None:
        clusters.append(current)
    return clusters


def seq_id_to_info(seq_id, strain_chrom_map):
    """Extract strain, chrom, detection from a sequence ID."""
    # nhmmer: STRAIN__chrIX_start_end_nhmmer
    if "__chr" in seq_id:
        parts = seq_id.split("__")
        strain = parts[0]
        chrom_part = parts[1]
        chrom = re.match(r'(chr[IVXLCDM]+)', chrom_part)
        chrom = chrom.group(1) if chrom else "?"
        det = "nhmmer" if seq_id.endswith("_nhmmer") else "BLAST"
    else:
        # plain strain name
        strain = seq_id
        chrom = strain_chrom_map.get(strain, "?")
        det = "BLAST"
    return strain, chrom, det


def build_cluster_df(gene, threshold_pct, strain_chrom_map):
    clstr_path = CDHIT_DIR / f"{gene}_cdhit_{threshold_pct}.clstr"
    if not clstr_path.exists():
        return pd.DataFrame()
    clusters = parse_clstr(clstr_path)
    rows = []
    for cl in clusters:
        members_detail = []
        chroms = set()
        species_set = set()
        detection_set = set()
        for seq_id in cl["members"]:
            strain, chrom, det = seq_id_to_info(seq_id, strain_chrom_map)
            spp = s2spp.get(strain, "?").replace("Saccharomyces ","S. ")
            chroms.add(chrom)
            species_set.add(spp)
            detection_set.add(det)
            members_detail.append(f"{strain}({chrom})")
        rows.append({
            "cluster": cl["id"],
            "n_members": len(cl["members"]),
            "representative": cl["representative"],
            "chromosomes": "; ".join(sorted(chroms)),
            "species": "; ".join(sorted(species_set)),
            "detection": "; ".join(sorted(detection_set)),
            "members": " | ".join(members_detail),
        })
    return pd.DataFrame(rows)


print("=" * 70)
print("CD-HIT CLUSTER ANALYSIS  (with chromosome mapping)")
print("=" * 70)

for gene in ("KHS1", "KHR1"):
    strain_chrom = load_strain_chrom(gene)
    print(f"\n{'─'*55}")
    print(f"  {gene}  (strain→chrom map: {len(strain_chrom)} entries)")
    print(f"{'─'*55}")
    for t in (90, 80, 70, 60, 50):
        df = build_cluster_df(gene, t, strain_chrom)
        if df.empty:
            continue
        out_tsv = CDHIT_DIR / f"{gene}_clusters_{t}pct.tsv"
        df.to_csv(out_tsv, sep="\t", index=False)
        print(f"\n  @ {t}%: {len(df)} clusters")
        for _, row in df.iterrows():
            print(f"    C{row['cluster']:02d} ({row['n_members']:2d} seqs) | "
                  f"chrom: {row['chromosomes']:<25} | {row['species'][:50]}")

print("\nAll cluster TSVs saved to:", CDHIT_DIR)
