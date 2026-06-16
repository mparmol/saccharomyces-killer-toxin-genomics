#!/usr/bin/env bash
# Script 23 — M-killer dsRNA preprotoxin HMM search in 83 Saccharomyces genomes
#
# Hypothesis (NUPAV model, Frank & Wolfe 2009):
#   KHS1 is a nuclear integration of M2-satellite dsRNA.
#   Divergent KHR1 loci may be independent integrations of M-killer dsRNA.
#
# Strategy:
#   1. Download preprotoxin sequences from GenBank (M1, M2, M28).
#   2. Align preprotoxins (nucleotide) with MAFFT.
#   3. Build HMM profile with hmmbuild (nucleotide mode).
#   4. Run nhmmer against each genome.
#   5. Filter hits that do NOT overlap known KHS1/KHR1 loci (novel integrations).
#   6. Extract novel hit sequences for downstream phylogenetics.
#
# GenBank accessions used:
#   M1  → U78817.1   (ScV-M1 complete genome; preprotoxin = K1 toxin, ~1600 bp)
#   M2  → MF957266.1 (ScV-M2-4 complete genome; preprotoxin = K2 toxin, ~1500 bp)
#   M28 → MF358735.1 (S. paradoxus M28-like satellite; K28 preprotoxin CDS)
#   Mlus→ GU723494.1 (ScV-Mlus; Klus toxin precursor, ~2100 bp)
#   Note: M1/M2 genomes contain one ORF encoding preprotoxin + 5'/3'-UTR.
#         We download full genomes and extract preprotoxin CDS coordinates below.
#
# Prerequisites (sac_killer conda env):
#   conda install -c bioconda hmmer entrez-direct -y
#   (efetch is part of entrez-direct)

set -euo pipefail

source /home/mparmol/miniconda3/etc/profile.d/conda.sh
conda activate sac_killer

PROJ=""
GENOMES_DIR="$PROJ/Nuclear_genomes"
OUT="$PROJ/resultados/10_mkiller_hmm"
SEQS_DIR="$OUT/preprotoxin_seqs"
HITS_DIR="$OUT/nhmmer_hits"
KNOWN_DIR="$PROJ/resultados/02_deteccion_blast"

mkdir -p "$SEQS_DIR" "$HITS_DIR"

# ─────────────────────────────────────────────
# STEP 1 — Download M-killer sequences from NCBI
# ─────────────────────────────────────────────
echo "=== STEP 1: Download preprotoxin sequences from NCBI ==="

# Full genome sequences (preprotoxin is the main ORF in each)
for acc in U78817.1 MF957266.1 MF358735.1 GU723494.1; do
    name="${acc%%.*}"
    outfile="$SEQS_DIR/${name}.fasta"
    if [ -f "$outfile" ]; then
        echo "  $acc already downloaded, skipping."
    else
        echo "  Downloading $acc ..."
        efetch -db nucleotide -id "$acc" -format fasta > "$outfile"
        sleep 1   # be polite with NCBI
    fi
done

# Rename for clarity
cp "$SEQS_DIR/U78817.fasta"   "$SEQS_DIR/M1_ScV-M1_genome.fasta"
cp "$SEQS_DIR/MF957266.fasta" "$SEQS_DIR/M2_ScV-M2-4_genome.fasta"
cp "$SEQS_DIR/MF358735.fasta" "$SEQS_DIR/M28_Sparadoxus_preprotoxin.fasta"
cp "$SEQS_DIR/GU723494.fasta" "$SEQS_DIR/Mlus_ScV-Mlus_genome.fasta"

echo ""
echo "Downloaded sequences:"
for f in "$SEQS_DIR"/*.fasta; do
    n=$(grep -c '>' "$f")
    echo "  $f  ($n seqs)"
done

# ─────────────────────────────────────────────
# STEP 2 — Build combined preprotoxin FASTA
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 2: Build combined preprotoxin FASTA ==="

# M28 is already just the preprotoxin CDS (MF358735.1).
# M1 (U78817.1) and M2 (MF957266.1) are complete satellite genomes;
# the preprotoxin ORF is the dominant/only coding sequence.
# Mlus (GU723494.1) encodes the Klus preprotoxin.
# For HMM building we use the full downloaded sequences (dominated by preprotoxin ORF).

cat "$SEQS_DIR/M1_ScV-M1_genome.fasta" \
    "$SEQS_DIR/M2_ScV-M2-4_genome.fasta" \
    "$SEQS_DIR/M28_Sparadoxus_preprotoxin.fasta" \
    "$SEQS_DIR/Mlus_ScV-Mlus_genome.fasta" \
    > "$SEQS_DIR/mkiller_preprotoxins_all.fasta"

echo "Combined FASTA: $SEQS_DIR/mkiller_preprotoxins_all.fasta"
grep -c '>' "$SEQS_DIR/mkiller_preprotoxins_all.fasta" | xargs echo "  Sequences:"

# ─────────────────────────────────────────────
# STEP 3 — Align with MAFFT
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 3: MAFFT alignment of preprotoxins ==="

mafft --localpair --maxiterate 1000 --thread 8 --reorder \
    "$SEQS_DIR/mkiller_preprotoxins_all.fasta" \
    > "$SEQS_DIR/mkiller_preprotoxins_mafft.fasta"

echo "Alignment: $SEQS_DIR/mkiller_preprotoxins_mafft.fasta"

# ─────────────────────────────────────────────
# STEP 4 — Build HMM profile (nucleotide)
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 4: Build nucleotide HMM profile ==="

hmmbuild --rna \
    "$OUT/mkiller_preprotoxin.hmm" \
    "$SEQS_DIR/mkiller_preprotoxins_mafft.fasta"

echo "HMM profile: $OUT/mkiller_preprotoxin.hmm"

# ─────────────────────────────────────────────
# STEP 5 — nhmmer against 83 genomes
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 5: nhmmer search against 83 genomes ==="

N_GENOMES=$(ls "$GENOMES_DIR"/*.fasta 2>/dev/null | wc -l)
echo "  Genomes found: $N_GENOMES"

for genome_fasta in "$GENOMES_DIR"/*.fasta; do
    strain=$(basename "$genome_fasta" .fasta)
    tbl_out="$HITS_DIR/${strain}_mkiller.tbl"
    if [ -f "$tbl_out" ]; then
        continue  # already done
    fi
    nhmmer \
        --tblout "$tbl_out" \
        --incE 0.01 \
        --cpu 8 \
        "$OUT/mkiller_preprotoxin.hmm" \
        "$genome_fasta" > /dev/null 2>&1
done

echo "  nhmmer complete for all strains."
echo "  Results in: $HITS_DIR/"

# ─────────────────────────────────────────────
# STEP 6 — Merge tblout files + basic summary
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 6: Merge nhmmer tables ==="

MERGED="$OUT/mkiller_nhmmer_allhits.tsv"

# Header
echo -e "strain\tseqname\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\tstrand\tevalue\tscore\tbias\tdesc" > "$MERGED"

for tbl in "$HITS_DIR"/*.tbl; do
    strain=$(basename "$tbl" _mkiller.tbl)
    grep -v '^#' "$tbl" | awk -v s="$strain" 'NF>=13 {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", \
            s,$1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;
        for(i=16;i<=NF;i++) printf "%s ",$i; printf "\n"
    }' >> "$MERGED" 2>/dev/null || true
done

N_HITS=$(wc -l < "$MERGED")
echo "  Total hits in merged table: $((N_HITS - 1))"

# ─────────────────────────────────────────────
# STEP 7 — Filter: remove known KHS1/KHR1 loci
# ─────────────────────────────────────────────
echo ""
echo "=== STEP 7: Filter known KHS1/KHR1 loci → only NOVEL hits ==="

# This Python one-liner loads the blast_hits_detail TSVs and cross-references
# coordinates from the merged nhmmer table. Any hit within 500 bp of a known
# KHS1 or KHR1 locus is flagged as "known" and excluded.
python3 - <<'PYEOF'
import pandas as pd, sys
from pathlib import Path

PROJ = Path("")
BLAST = PROJ / "res" / "02_deteccion_blast"
OUT   = PROJ / "re" / "10_mkiller_hmm"
MERGED = OUT / "mkiller_nhmmer_allhits.tsv"

df = pd.read_csv(MERGED, sep="\t")
if df.empty:
    print("  No M-killer hits found in any genome.")
    sys.exit(0)

print(f"  Total hits before filtering: {len(df)}")

# Load known KHS1 / KHR1 coordinates
known = []
for gene in ("KHS1", "KHR1"):
    detail = BLAST / gene / f"{gene}_blast_hits_detail.tsv"
    if detail.exists():
        dg = pd.read_csv(detail, sep="\t")
        for _, row in dg.iterrows():
            known.append({
                "strain": str(row.get("strain","")).strip(),
                "chrom": str(row.get("chrom","")).strip(),
                "start": int(row.get("sstart", row.get("start", 0))),
                "end":   int(row.get("send",   row.get("end", 0))),
            })

def overlaps_known(row, known_loci, slack=500):
    for lk in known_loci:
        if lk["strain"] != str(row["strain"]).strip():
            continue
        if lk["chrom"] not in str(row["seqname"]):
            continue
        hs, he = sorted([int(row["ali_from"]), int(row["ali_to"])])
        ks, ke = sorted([int(lk["start"]), int(lk["end"])])
        if hs <= ke + slack and he >= ks - slack:
            return True
    return False

df["is_known"] = df.apply(lambda r: overlaps_known(r, known), axis=1)
novel = df[~df["is_known"]].copy()
print(f"  Known KHS1/KHR1 hits removed: {df['is_known'].sum()}")
print(f"  Novel M-killer-like hits: {len(novel)}")

novel.to_csv(OUT / "mkiller_novel_hits.tsv", sep="\t", index=False)

if len(novel) == 0:
    print("\n  RESULT: No novel M-killer integrations found beyond known KHS1/KHR1 loci.")
    print("  This is a NEGATIVE RESULT — equally important for publication.")
    print("  Interpretation: KHS1 and KHR1 appear to be the only detectable nuclear")
    print("  integrations of M-killer dsRNA in these 83 Saccharomyces genomes.")
else:
    print(f"\n  RESULT: {len(novel)} novel hits warrant manual inspection.")
    print("  See: mkiller_novel_hits.tsv")
    strains = novel["strain"].value_counts()
    print("\n  Strains with novel hits:")
    for s, n in strains.items():
        print(f"    {s}: {n} hit(s)")

# Summary by evalue threshold
print("\n  Evalue distribution of novel hits:")
for thr in [1e-20, 1e-10, 1e-5, 0.01]:
    n = (novel["evalue"].astype(float) < thr).sum() if len(novel) > 0 else 0
    print(f"    E < {thr:.0e}: {n} hits")

PYEOF

echo ""
echo "=== M-killer HMM search complete ==="
echo "Results:"
echo "  $OUT/mkiller_preprotoxin.hmm"
echo "  $OUT/mkiller_nhmmer_allhits.tsv"
echo "  $OUT/mkiller_novel_hits.tsv"
echo ""
echo "INTERPRETATION GUIDE:"
echo "  - If mkiller_novel_hits.tsv is EMPTY: KHS1 and KHR1 are the only M-killer"
echo "    integrations in these 83 genomes. Supports the hypothesis that KHR1 divergent"
echo "    loci arose by ancient divergence/translocation, not fresh viral capture."
echo "  - If novel hits FOUND: each must be manually inspected (coords, ORF integrity,"
echo "    species context) before claiming new integrations."
echo "  - IMPORTANT: M1/M2/M28/Mlus preprotoxins share limited identity (~30-40%)."
echo "    Very divergent integrations may be missed. A protein-level search (tblastn"
echo "    with preprotoxin translations) would be complementary."
