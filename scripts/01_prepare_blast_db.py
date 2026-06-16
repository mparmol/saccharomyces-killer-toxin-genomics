#!/usr/bin/env python3
"""
Step 1: Merge all 83 nuclear genomes into a single FASTA with strain-prefixed headers,
then build a BLAST nucleotide database.
"""

import os
import subprocess
import sys
from pathlib import Path

PROJECT_DIR = Path("")
GENOMES_DIR = PROJECT_DIR / "Nuclear genomes"
RESULTS_DIR = PROJECT_DIR / "results"
BLAST_DB_DIR = RESULTS_DIR / "blast_db"

MERGED_FASTA = BLAST_DB_DIR / "all_genomes_merged.fasta"
BLAST_DB     = BLAST_DB_DIR / "all_genomes"


def merge_genomes():
    genome_files = sorted(GENOMES_DIR.glob("*.fasta"))
    print(f"Found {len(genome_files)} genome assemblies")

    total_seqs = 0
    with open(MERGED_FASTA, "w") as outf:
        for gf in genome_files:
            strain = gf.stem
            with open(gf) as inf:
                for line in inf:
                    if line.startswith(">"):
                        chrom = line[1:].strip().split()[0]
                        outf.write(f">{strain}__{chrom}\n")
                        total_seqs += 1
                    else:
                        outf.write(line)

    print(f"Merged {total_seqs} sequences from {len(genome_files)} genomes → {MERGED_FASTA}")


def build_blast_db():
    cmd = [
        "makeblastdb",
        "-in",  str(MERGED_FASTA),
        "-dbtype", "nucl",
        "-out", str(BLAST_DB),
        "-parse_seqids",
        "-title", "All_Saccharomyces_genomes",
    ]
    print("Building BLAST nucleotide database...")
    subprocess.run(cmd, check=True)
    print(f"BLAST DB ready: {BLAST_DB}")


if __name__ == "__main__":
    merge_genomes()
    build_blast_db()
