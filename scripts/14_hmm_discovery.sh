#!/usr/bin/env bash

set -euo pipefail

PROJECT=""
GENOME_DIR="$PROJECT/Nuclear genomes"
ALN_DIR="$PROJECT/results/alignments"
DISC_DIR="$PROJECT/results_2/discovery"
COMB_DIR="$DISC_DIR/combined"
PFAM_DB="$PROJECT/databases/pfam/Pfam-A.hmm"
EGGNOG_DB="$PROJECT/databases/eggnog"
THREADS="${1:-8}"

CONDA_BASE="${CONDA_BASE:-}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate sac_killer


mkdir -p "$DISC_DIR/profiles" "$DISC_DIR/hits" "$DISC_DIR/annotation" "$COMB_DIR"

echo ""

hmmbuild \
    "$DISC_DIR/profiles/KHS1.hmm" \
    "$ALN_DIR/KHS1/KHS1_primary_mafft_trim.fasta"
echo "  KHS1.hmm construido"

hmmbuild \
    "$DISC_DIR/profiles/KHR1.hmm" \
    "$ALN_DIR/KHR1/KHR1_mafft_gt50.fasta"
echo "  KHR1.hmm construido"

echo ""

for GENE in KHS1 KHR1; do
    HITS_DIR="$DISC_DIR/hits/$GENE"
    mkdir -p "$HITS_DIR"
    echo "  Gene: $GENE"

    run_nhmmer() {
        local fasta="$1"
        local gene="$2"
        local hits_dir="$3"
        local strain=$(basename "$fasta" .fasta)
        local out="$hits_dir/${strain}.tbl"
        if [ -f "$out" ]; then
            return  
        fi
        nhmmer \
            --tblout "$out" \
            --noali \
            --incE 0.01 \
            --cpu 1 \
            "$PROJECT/results_2/discovery/profiles/${gene}.hmm" \
            "$fasta" 2>/dev/null
    }
    export -f run_nhmmer
    export PROJECT

    find "$GENOME_DIR" -name "*.fasta" | \
        xargs -P "$THREADS" -I{} bash -c "run_nhmmer '{}' '$GENE' '$HITS_DIR'"

done

echo ""

python3 "$PROJECT/scripts/14b_process_hits.py" \
    --project_dir "$PROJECT" \
    --threads "$THREADS"

echo ""

for GENE in KHS1 KHR1; do
    PROT_ALL="$DISC_DIR/hits/${GENE}_all_proteins.fasta"
    PROT_NEW="$DISC_DIR/hits/${GENE}_new_candidates_protein.fasta"

    if [ ! -f "$PROT_ALL" ]; then
        continue
    fi

    hmmscan \
        --domtblout "$DISC_DIR/annotation/pfam_${GENE}_domtbl.tsv" \
        --tblout     "$DISC_DIR/annotation/pfam_${GENE}_tbl.tsv" \
        -E 1e-3 \
        --cpu "$THREADS" \
        "$PFAM_DB" \
        "$PROT_ALL" \
        > "$DISC_DIR/annotation/pfam_${GENE}.log" 2>&1

    if [ -f "$PROT_NEW" ] && [ -s "$PROT_NEW" ]; then
        echo "  EggNOG — $GENE candidatos nuevos..."
        emapper.py \
            -i "$PROT_NEW" \
            --data_dir "$EGGNOG_DB" \
            -o "${GENE}_new_candidates" \
            --output_dir "$DISC_DIR/annotation" \
            --cpu "$THREADS" \
            --override 2>/dev/null
    else
    fi
done


for GENE in KHS1 KHR1; do
    NEW_SEQS="$DISC_DIR/hits/${GENE}_new_candidates_nucl.fasta"
    if [ ! -f "$NEW_SEQS" ] || [ ! -s "$NEW_SEQS" ]; then
        continue
    fi

    case $GENE in
        KHS1) BASE_ALN="$PROJECT/results/alignments/KHS1/KHS1_primary_mafft_trim.fasta" ;;
        KHR1) BASE_ALN="$PROJECT/results/alignments/KHR1/KHR1_mafft_gt50.fasta" ;;
    esac

    mafft \
        --add "$NEW_SEQS" \
        --keeplength \
        --reorder \
        --thread "$THREADS" \
        "$BASE_ALN" \
        > "$COMB_DIR/${GENE}_extended_mafft.fasta" 2>/dev/null

    trimal \
        -in  "$COMB_DIR/${GENE}_extended_mafft.fasta" \
        -out "$COMB_DIR/${GENE}_extended_trim.fasta" \
        -gt 0.5

    iqtree \
        -s "$COMB_DIR/${GENE}_extended_trim.fasta" \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "$THREADS" \
        --redo \
        --prefix "$COMB_DIR/${GENE}_extended_iqtree" \
        2>/dev/null

done

