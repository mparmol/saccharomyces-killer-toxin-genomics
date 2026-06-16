#!/usr/bin/env python3
"""
Step 2: Run BLASTn and tBLASTx for KHS1 and KHR1 against the merged genome database.

Strategy:
  - BLASTn: fast, high identity — catches close homologs
  - tBLASTx: translated 6-frame search — catches divergent homologs
  - Both are run for completeness; sequences from both are merged and deduplicated
"""

import subprocess
from pathlib import Path

PROJECT_DIR  = Path("")
RESULTS_DIR  = PROJECT_DIR / "results"
BLAST_DB     = RESULTS_DIR / "blast_db" / "all_genomes"
BLAST_OUT    = RESULTS_DIR / "blast_results"

QUERIES = {
    "KHS1": PROJECT_DIR / "KHS1.fas",
    "KHR1": PROJECT_DIR / "SK1_KHR1.fas",
}

THREADS = 8

# outfmt 6 columns (custom)
OUTFMT = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"


def run_blastn(gene, query, out_dir, evalue="1e-5", perc_id=60):
    outfile = out_dir / f"{gene}_blastn.tsv"
    cmd = [
        "blastn",
        "-query",       str(query),
        "-db",          str(BLAST_DB),
        "-out",         str(outfile),
        "-outfmt",      OUTFMT,
        "-evalue",      evalue,
        "-perc_identity", str(perc_id),
        "-num_threads", str(THREADS),
        "-word_size",   "11",
    ]
    print(f"  Running BLASTn for {gene}...")
    subprocess.run(cmd, check=True)
    print(f"  → {outfile}")
    return outfile


def run_tblastx(gene, query, out_dir, evalue="1e-5"):
    outfile = out_dir / f"{gene}_tblastx.tsv"
    cmd = [
        "tblastx",
        "-query",       str(query),
        "-db",          str(BLAST_DB),
        "-out",         str(outfile),
        "-outfmt",      OUTFMT,
        "-evalue",      evalue,
        "-num_threads", str(THREADS),
    ]
    print(f"  Running tBLASTx for {gene}...")
    subprocess.run(cmd, check=True)
    print(f"  → {outfile}")
    return outfile


if __name__ == "__main__":
    for gene, query in QUERIES.items():
        out_dir = BLAST_OUT / gene
        print(f"\n=== {gene} ===")
        run_blastn(gene, query, out_dir)
        run_tblastx(gene, query, out_dir)

    print("\nBLAST searches complete.")
