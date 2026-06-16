
import os
import re
import pandas as pd

PROJECT    = ""
GENOME_DIR = os.path.join(PROJECT, "Nuclear genomes")
EXCEL      = os.path.join(PROJECT, "Table Strain Information killer.xlsx")
OUT_DIR    = os.path.join(PROJECT, "results")
OUT_TSV    = os.path.join(OUT_DIR, "genome_metadata.tsv")
OUT_XLSX   = os.path.join(OUT_DIR, "genome_metadata.xlsx")

# ── 1. Leer Excel ────────────────────────────────────────────────────────────
df_meta = pd.read_excel(EXCEL, header=2)
df_meta.columns = [
    "Strain", "Species", "NCBI_Biosample", "Synonym",
    "Population", "Reference", "Nuclear_Accession", "KHR1", "KHS1"
]
df_meta = df_meta.dropna(subset=["Strain"]).copy()
df_meta["Strain"] = df_meta["Strain"].astype(str).str.strip()

def normalize(s):
    return re.sub(r'[-.\s]', '', s).upper()

STRAIN_CORRECTIONS = {
    "UFRJ50816T": "UFRJ50816",   
}
df_meta["strain_key"] = df_meta["Strain"].apply(
    lambda s: normalize(STRAIN_CORRECTIONS.get(s, s))
)

def parse_fasta(path):
    """Devuelve lista de (header, seq_str) desde un FASTA."""
    seqs = []
    header = None
    buf = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if header is not None:
        seqs.append((header, "".join(buf)))
    return seqs

def gc_content(seqs):
    total_gc = sum(s.count('G') + s.count('C') for _, s in seqs)
    total_acgt = sum(sum(s.count(b) for b in 'ACGT') for _, s in seqs)
    return round(100 * total_gc / total_acgt, 2) if total_acgt else 0.0

def n50(lengths):
    sl = sorted(lengths, reverse=True)
    half = sum(sl) / 2
    cum = 0
    for l in sl:
        cum += l
        if cum >= half:
            return l
    return 0

rows = []
fasta_files = sorted(f for f in os.listdir(GENOME_DIR) if f.endswith(".fasta"))
print(f"Procesando {len(fasta_files)} genomas...")

for fname in fasta_files:
    strain_file = fname.replace(".fasta", "")
    fpath = os.path.join(GENOME_DIR, fname)

    seqs    = parse_fasta(fpath)
    lengths = [len(s) for _, s in seqs]
    total   = sum(lengths)
    gc      = gc_content(seqs)
    n50v    = n50(lengths)

    rows.append({
        "Strain_file":     strain_file,
        "strain_key":      normalize(strain_file),
        "N_chromosomes":   len(seqs),
        "Genome_size_bp":  total,
        "Genome_size_Mb":  round(total / 1e6, 3),
        "GC_percent":      gc,
        "N50_bp":          n50v,
        "Largest_chr_bp":  max(lengths),
        "Smallest_chr_bp": min(lengths),
        "Chromosomes":     ";".join(h for h, _ in seqs),
    })
    print(f"  {strain_file}: {len(seqs)} chr, {total/1e6:.3f} Mb, GC={gc}%")

df_fasta = pd.DataFrame(rows)

df = pd.merge(df_meta, df_fasta, on="strain_key", how="outer")

no_fasta = df[df["Strain_file"].isna()]["Strain"].tolist()
no_meta  = df[df["Strain"].isna()]["Strain_file"].tolist()
if no_fasta:
    print(f"\n: {no_fasta}")
if no_meta:
    print(f": {no_meta}")

cols = [
    "Strain", "Species", "Population", "NCBI_Biosample", "Synonym",
    "Reference", "Nuclear_Accession",
    "N_chromosomes", "Genome_size_bp", "Genome_size_Mb",
    "GC_percent", "N50_bp", "Largest_chr_bp", "Smallest_chr_bp",
    "KHR1", "KHS1", "Chromosomes"
]
cols = [c for c in cols if c in df.columns]
df = df[cols].sort_values("Strain").reset_index(drop=True)

os.makedirs(OUT_DIR, exist_ok=True)
df.to_csv(OUT_TSV, sep="\t", index=False)

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    df.to_excel(writer, index=False, sheet_name="Genome_metadata")
    ws = writer.sheets["Genome_metadata"]
    for col in ws.columns:
        max_w = max((len(str(cell.value)) if cell.value else 0 for cell in col), default=10)
        ws.column_dimensions[col[0].column_letter].width = min(max_w + 2, 55)

print(f"\n: {OUT_TSV}")
print(f": {OUT_XLSX}")
print(f": {len(df)}")

print("\n")
summary = df.groupby("Species").agg(
    N_strains      = ("Strain", "count"),
    Mb_mean        = ("Genome_size_Mb", "mean"),
    Mb_min         = ("Genome_size_Mb", "min"),
    Mb_max         = ("Genome_size_Mb", "max"),
    GC_mean        = ("GC_percent", "mean"),
    N_chr_median   = ("N_chromosomes", "median"),
    N50_mean_kb    = ("N50_bp", lambda x: round(x.mean() / 1e3, 1)),
).round(3)
print(summary.to_string())
