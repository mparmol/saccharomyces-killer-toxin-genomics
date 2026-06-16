#!/usr/bin/env bash

set -euo pipefail

PROJECT=""
GENOME_DIR="$PROJECT/Nuclear genomes"
BUSCO_OUT="$PROJECT/results_2/phylogenomics/busco"
BUSCO_DB="$PROJECT/databases/busco"
OUT_DIR="$PROJECT/results_2/phylogenomics"
THREADS="${1:-8}"

CONDA_BASE="${CONDA_BASE:}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate sac_killer


echo ""
echo "[1/4] ..."
mkdir -p "$BUSCO_OUT"

run_busco() {
    local fasta="$1"
    local strain=$(basename "$fasta" .fasta)
    local out="$BUSCO_OUT/$strain"
    if [ -f "$out/run_saccharomycetes_odb10/full_table.tsv" ]; then
        echo "  [OK] $strain — ya procesado"
        return
    fi
    echo "  $strain..."
    busco \
        -m genome \
        -i "$fasta" \
        -o "$strain" \
        --out_path "$BUSCO_OUT" \
        -l saccharomycetes_odb10 \
        --download_path "$BUSCO_DB" \
        --cpu $((THREADS / 4 < 1 ? 1 : THREADS / 4)) \
        --offline \
        --skip_bbtools \
        -q 2>/dev/null
}
export -f run_busco
export BUSCO_OUT BUSCO_DB THREADS

find "$GENOME_DIR" -name "*.fasta" | \
    xargs -P 4 -I{} bash -c 'run_busco "$@"' _ {}


echo ""

python3 "$PROJECT/scripts/13b_build_supermatrix.py" \
    --busco_dir "$BUSCO_OUT" \
    --genome_dir "$GENOME_DIR" \
    --out_dir "$OUT_DIR" \
    --min_occupancy 0.96 \
    --threads "$THREADS"

SUPERMATRIX="$OUT_DIR/supermatrix/supermatrix.fasta"
PARTITIONS="$OUT_DIR/supermatrix/partitions.nex"
PREFIX="$OUT_DIR/supermatrix/phylogenomics_iqtree"

if [ ! -f "$SUPERMATRIX" ]; then
    exit 1
fi

N_GENES=$(grep -c "^>" "$SUPERMATRIX" | head -1 || true)
echo ""

iqtree \
    -s "$SUPERMATRIX" \
    -p "$PARTITIONS" \
    -m MFP+MERGE \
    -B 1000 \
    --alrt 1000 \
    -T "$THREADS" \
    --redo \
    --prefix "$PREFIX"


echo ""
python3 - <<'PYEOF'
import os, glob, re
from pathlib import Path

busco_dir = Path(os.environ["BUSCO_OUT"])
rows = []
for summary in sorted(busco_dir.glob("*/short_summary.specific.saccharomycetes_odb10.*.txt")):
    strain = summary.parts[-3]
    text = summary.read_text()
    m = re.search(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%", text)
    if m:
        rows.append(f"{strain}\tC:{m.group(1)}%  S:{m.group(2)}%  D:{m.group(3)}%  F:{m.group(4)}%  M:{m.group(5)}%")

summary_out = busco_dir.parent / "busco_summary.tsv"
summary_out.write_text("Strain\tBUSCO_completeness\n" + "\n".join(rows))
print(f"  Resumen BUSCO: {summary_out}")
PYEOF