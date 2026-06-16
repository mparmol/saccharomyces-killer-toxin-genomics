#!/usr/bin/env bash
# Step 4: Multiple sequence alignment (MAFFT) + trimming (trimAl) + phylogenetic tree (IQ-TREE2)
#
# For each gene (KHS1, KHR1):
#   1. MAFFT --auto alignment
#   2. trimAl -automated1 trimming (removes gappy columns)
#   3. IQ-TREE2: ModelFinder + UFBoot 1000 + SH-aLRT 1000
#

set -euo pipefail

PROJECT=""
SEQ_DIR="$PROJECT/results/sequences"
ALN_DIR="$PROJECT/results/alignments"
TREE_DIR="$PROJECT/results/trees"
THREADS=8

for GENE in KHS1 KHR1; do
    echo "=============================="
    echo " Processing: $GENE"
    echo "=============================="

    INPUT="$SEQ_DIR/$GENE/${GENE}_hits.fasta"
    ALN_RAW="$ALN_DIR/$GENE/${GENE}_mafft_raw.fasta"
    ALN_TRIM="$ALN_DIR/$GENE/${GENE}_mafft_trim.fasta"
    TREE_OUT="$TREE_DIR/$GENE"

    # Skip if input doesn't exist or is empty
    if [[ ! -s "$INPUT" ]]; then
        echo "  WARN: No sequences found for $GENE, skipping."
        continue
    fi

    N_SEQS=$(grep -c "^>" "$INPUT")
    echo "  Input sequences: $N_SEQS"

    # --- 1. MAFFT alignment ---
    echo "  Running MAFFT..."
    if [[ $N_SEQS -le 200 ]]; then
        # L-INS-i: most accurate for < 200 sequences
        mafft --localpair --maxiterate 1000 --thread $THREADS \
              --reorder --adjustdirection \
              "$INPUT" > "$ALN_RAW" 2>"$ALN_DIR/$GENE/${GENE}_mafft.log"
    else
        mafft --auto --thread $THREADS --reorder --adjustdirection \
              "$INPUT" > "$ALN_RAW" 2>"$ALN_DIR/$GENE/${GENE}_mafft.log"
    fi
    echo "  MAFFT done → $ALN_RAW"

    # --- 2. trimAl trimming ---
    echo "  Trimming with trimAl..."
    trimal \
        -in  "$ALN_RAW" \
        -out "$ALN_TRIM" \
        -automated1 \
        -htmlout "$ALN_DIR/$GENE/${GENE}_trimal.html" \
        -colnumbering 2>"$ALN_DIR/$GENE/${GENE}_trimal.log"
    echo "  trimAl done → $ALN_TRIM"

    # --- 3. IQ-TREE2 ---
    echo "  Running IQ-TREE2..."
    mkdir -p "$TREE_OUT"
    iqtree \
        -s       "$ALN_TRIM" \
        --prefix "$TREE_OUT/${GENE}_iqtree" \
        -m       MFP \
        -B       1000 \
        --alrt   1000 \
        -T       $THREADS \
        --redo \
        2>&1 | tee "$TREE_OUT/${GENE}_iqtree.log"

    echo "  IQ-TREE2 done. Best tree: $TREE_OUT/${GENE}_iqtree.treefile"
    echo ""
done

echo "All done. Trees in $TREE_DIR"
