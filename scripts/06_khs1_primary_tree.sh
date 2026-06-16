#!/usr/bin/env bash
# Build a separate tree for KHS1 using only PRIMARY copy sequences
# (excludes _R_ reversed sequences which form a paralog group).
# This gives a cleaner per-strain phylogeny.

set -euo pipefail

PROJECT=""
SEQ_DIR="$PROJECT/results/sequences/KHS1"
ALN_DIR="$PROJECT/results/alignments/KHS1"
TREE_DIR="$PROJECT/results/trees/KHS1"
THREADS=8

echo "=== KHS1 PRIMARY-ONLY TREE ==="

# Filter out _R_ sequences and keep only copy1 from strains with 2 copies
# (copy1 = the longer/better hit, copy2 = the paralog that got reverse-complemented by MAFFT)
python3 - <<'EOF'
from Bio import SeqIO
from pathlib import Path

raw_fasta  = Path("KHS1_hits.fasta")
trim_aln   = Path("KHS1/KHS1_mafft_trim.fasta")
out_fasta  = Path("KHS1/KHS1_primary_hits.fasta")

# Read all raw hit sequences; exclude _copy2 (which are the paralogous shorter hits)
# Primary = sequences without _copy suffix or with _copy1 but NOT _copy2
recs = []
for rec in SeqIO.parse(raw_fasta, "fasta"):
    sid = rec.id
    # Exclude copy2 and higher (paralog group)
    if "_copy2" in sid or "_copy3" in sid:
        continue
    # Rename _copy1 to just the strain name
    rec.id = sid.replace("_copy1", "")
    rec.description = ""
    recs.append(rec)

SeqIO.write(recs, out_fasta, "fasta")
print(f"Primary KHS1 sequences: {len(recs)} -> {out_fasta}")
EOF

INPUT="$SEQ_DIR/KHS1_primary_hits.fasta"
ALN_RAW="$ALN_DIR/KHS1_primary_mafft_raw.fasta"
ALN_TRIM="$ALN_DIR/KHS1_primary_mafft_trim.fasta"
TREE_OUT="$TREE_DIR"
N_SEQS=$(grep -c "^>" "$INPUT")
echo "  Input: $N_SEQS sequences"

echo "  Running MAFFT (L-INS-i)..."
mafft --localpair --maxiterate 1000 --thread $THREADS --reorder \
      "$INPUT" > "$ALN_RAW" 2>"$ALN_DIR/KHS1_primary_mafft.log"

echo "  Trimming with trimAl..."
trimal -in "$ALN_RAW" -out "$ALN_TRIM" -automated1 \
       -htmlout "$ALN_DIR/KHS1_primary_trimal.html" \
       -colnumbering 2>"$ALN_DIR/KHS1_primary_trimal.log"

echo "  Running IQ-TREE2..."
iqtree \
    -s       "$ALN_TRIM" \
    --prefix "$TREE_OUT/KHS1_primary_iqtree" \
    -m       MFP \
    -B       1000 \
    --alrt   1000 \
    -T       $THREADS \
    --redo \
    2>&1 | tee "$TREE_OUT/KHS1_primary_iqtree.log"

echo "Done. Primary KHS1 tree: $TREE_OUT/KHS1_primary_iqtree.treefile"
