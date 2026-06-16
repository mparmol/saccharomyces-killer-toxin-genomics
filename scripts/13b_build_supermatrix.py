#!/usr/bin/env python3


import argparse
import os
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--busco_dir",      required=True, type=Path)
    p.add_argument("--genome_dir",     required=True, type=Path)
    p.add_argument("--out_dir",        required=True, type=Path)
    p.add_argument("--min_occupancy",  type=float, default=0.96)
    p.add_argument("--threads",        type=int, default=8)
    return p.parse_args()


def parse_full_table(tsv_path):
    """"""
    genes = {}
    for line in Path(tsv_path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        busco_id, status = parts[0], parts[1]
        seq_file = parts[2] if len(parts) > 2 else ""
        genes[busco_id] = (status, seq_file)
    return genes


def run(cmd, **kw):
    subprocess.run(cmd, shell=True, check=True, **kw)


def main():
    args = parse_args()
    busco_dir    = args.busco_dir
    out_dir      = args.out_dir
    orth_dir     = out_dir / "orthologs"
    aln_dir      = out_dir / "alignments"
    matrix_dir   = out_dir / "supermatrix"
    for d in [orth_dir, aln_dir, matrix_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("")
    strain_genes = {}   # {strain: {busco_id: fasta_path}}
    for strain_dir in sorted(busco_dir.iterdir()):
        if not strain_dir.is_dir():
            continue
        tables = list(strain_dir.glob("run_saccharomycetes_odb10/full_table.tsv"))
        if not tables:
            print(f"  [WARN] Sin full_table.tsv en {strain_dir.name}")
            continue
        strain = strain_dir.name
        gene_info = parse_full_table(tables[0])
        single_copy = {}
        for gid, (status, seq_file) in gene_info.items():
            if status == "Complete":
                # BUSCO 6 outputs .faa (protein) — no .fna nucleotide files
                faa_candidates = list(strain_dir.glob(
                    f"run_saccharomycetes_odb10/busco_sequences/single_copy_busco_sequences/{gid}.faa"
                ))
                if faa_candidates:
                    single_copy[gid] = faa_candidates[0]
        strain_genes[strain] = single_copy
        print(f"  {strain}: {len(single_copy)} genes single-copy")

    if not strain_genes:
        print("ERROR")
        sys.exit(1)

    n_strains = len(strain_genes)
    min_strains = int(args.min_occupancy * n_strains)
    print(f"\n : {n_strains}")
    print(f"  : {min_strains}/{n_strains} ({args.min_occupancy*100:.0f}%)")

    gene_counts = defaultdict(int)
    for sg in strain_genes.values():
        for gid in sg:
            gene_counts[gid] += 1

    shared_genes = [g for g, c in gene_counts.items() if c >= min_strains]
    shared_genes.sort()
    print(f"\n  Genes >= {min_strains}: {len(shared_genes)}")

    if len(shared_genes) < 10:
        print("")
        sys.exit(1)

    print(f"\n {len(shared_genes)} ...")
    strains_ordered = sorted(strain_genes.keys())
    aligned_files = []

    for i, gid in enumerate(shared_genes, 1):
        orth_fasta = orth_dir / f"{gid}.fna"
        aln_fasta  = aln_dir  / f"{gid}_trim.fasta"

        with open(orth_fasta, "w") as fout:
            for strain in strains_ordered:
                sg = strain_genes.get(strain, {})
                if gid in sg:
                    seq_path = sg[gid]
                    seqs = open(seq_path).read().strip().split("\n")
                    fout.write(f">{strain}\n")
                    fout.write("\n".join(l for l in seqs if not l.startswith(">")) + "\n")

        # MAFFT + trimAl
        if not aln_fasta.exists():
            try:
                run(f"mafft --auto --thread {args.threads} --quiet "
                    f"'{orth_fasta}' > '{aln_fasta}.tmp' 2>/dev/null")
                run(f"trimal -in '{aln_fasta}.tmp' -out '{aln_fasta}' "
                    f"-automated1 2>/dev/null")
                Path(f"{aln_fasta}.tmp").unlink(missing_ok=True)
            except subprocess.CalledProcessError:
                print(f"  [WARN]  {gid}")
                continue

        if aln_fasta.stat().st_size > 100:
            aligned_files.append((gid, aln_fasta))

        if i % 50 == 0:
            print(f"  {i}/{len(shared_genes)} ...")

    print(f"  : {len(aligned_files)}")

    print(f"\n")

    def read_fasta(path):
        seqs = {}
        header = None
        buf = []
        for line in Path(path).read_text().splitlines():
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if header:
            seqs[header] = "".join(buf)
        return seqs

    gene_lengths = {}
    gene_seqs = {}
    for gid, aln_f in aligned_files:
        seqs = read_fasta(aln_f)
        if not seqs:
            continue
        aln_len = len(next(iter(seqs.values())))
        gene_lengths[gid] = aln_len
        gene_seqs[gid] = seqs

    supermatrix = {}
    for strain in strains_ordered:
        row = []
        for gid, aln_f in aligned_files:
            if gid not in gene_seqs:
                continue
            seqs = gene_seqs[gid]
            aln_len = gene_lengths[gid]
            seq = seqs.get(strain, "-" * aln_len)
            # Enforce uniform length (truncate or pad with gaps)
            if len(seq) > aln_len:
                seq = seq[:aln_len]
            elif len(seq) < aln_len:
                seq = seq + "-" * (aln_len - len(seq))
            row.append(seq)
        supermatrix[strain] = "".join(row)

    super_fasta = matrix_dir / "supermatrix.fasta"
    with open(super_fasta, "w") as f:
        for strain, seq in supermatrix.items():
            f.write(f">{strain}\n{seq}\n")
    total_len = sum(gene_lengths[g] for g in gene_lengths)
    print(f"  Supermatri: {len(supermatrix)} × {total_len} bp ({len(gene_lengths)} genes)")

    partition_nex = matrix_dir / "partitions.nex"
    pos = 1
    with open(partition_nex, "w") as f:
        f.write("#nexus\nbegin sets;\n")
        for gid, _ in aligned_files:
            if gid not in gene_lengths:
                continue
            l = gene_lengths[gid]
            f.write(f"  charset {gid} = {pos}-{pos + l - 1};\n")
            pos += l
        f.write("end;\n")
    print(f"  Partition file: {partition_nex}")
    print(f"  {len(gene_lengths)}")

    print(f"\n:")
    print(f"  iqtree2 -s {super_fasta} -p {partition_nex} -m MFP+MERGE -B 1000 -T {args.threads}")


if __name__ == "__main__":
    main()
