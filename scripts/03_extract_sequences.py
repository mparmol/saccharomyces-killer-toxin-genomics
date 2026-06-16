#!/usr/bin/env python3
"""
Step 3: Parse BLAST results and extract hit sequences from genome assemblies.

Strategy:
  - BLASTn hits: primary source (long, high-identity, cover most of CDS).
    Extract the hit region + FLANK bp on each side.
  - tBLASTx hits: used only for strains absent from BLASTn.
    tBLASTx produces many short overlapping hits (6-frame translation);
    hits on the same chromosome within CLUSTER_DIST bp are merged into one locus,
    then the merged span is extracted.
  - Per strain: keep best hit (highest bitscore). If multiple copies
    (bitscore >= PARALOG_FRAC × best), keep them all (numbered _copy1, _copy2 …).
"""

from pathlib import Path
from collections import defaultdict
import csv
from Bio import SeqIO

PROJECT_DIR = Path("")
GENOMES_DIR = PROJECT_DIR / "Nuclear genomes"
RESULTS_DIR = PROJECT_DIR / "results"
BLAST_DIR   = RESULTS_DIR / "blast_results"
SEQ_DIR     = RESULTS_DIR / "sequences"

GENES = ["KHS1", "KHR1"]

FLANK        = 150   # bp added on each side of extracted region
PARALOG_FRAC = 0.80  # keep copies with bitscore >= this fraction of best
CLUSTER_DIST = 500   # bp: merge tBLASTx hits closer than this
MIN_TBLASTX_BS     = 50   # ignore very weak tBLASTx hits
MIN_TBLASTX_PIDENT = 60   # require ≥60% aa identity for tBLASTx (prevents chrV false positives)


def parse_blast(tsv_file):
    cols = ["qseqid","sseqid","pident","length","qlen","slen",
            "qstart","qend","sstart","send","evalue","bitscore"]
    hits = []
    with open(tsv_file) as fh:
        for row in csv.DictReader(fh, fieldnames=cols, delimiter="\t"):
            for k in ["length","qlen","slen","qstart","qend","sstart","send"]:
                row[k] = int(row[k])
            row["pident"]   = float(row["pident"])
            row["evalue"]   = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            parts           = row["sseqid"].split("__", 1)
            row["strain"]   = parts[0]
            row["chrom"]    = parts[1] if len(parts) > 1 else parts[0]
            hits.append(row)
    return hits


def cluster_hits(hits, max_gap=CLUSTER_DIST):
    """
    Merge hits on the same (strain, chrom, strand) within max_gap bp.
    Returns list of merged loci dicts with: strain, chrom, start, end,
    strand ('+'/'-'), bitscore_sum, n_hits, pident_max, evalue_min.
    """
    # Normalise coordinates: start < end, record strand
    normed = []
    for h in hits:
        s, e = h["sstart"], h["send"]
        strand = "+" if s <= e else "-"
        normed.append({**h, "strand": strand,
                       "norm_start": min(s, e),
                       "norm_end":   max(s, e)})

    # Group by (strain, chrom, strand)
    groups = defaultdict(list)
    for h in normed:
        groups[(h["strain"], h["chrom"], h["strand"])].append(h)

    loci = []
    for (strain, chrom, strand), grp in groups.items():
        grp.sort(key=lambda x: x["norm_start"])
        cur_start = grp[0]["norm_start"]
        cur_end   = grp[0]["norm_end"]
        cur_hits  = [grp[0]]
        for h in grp[1:]:
            if h["norm_start"] - cur_end <= max_gap:
                cur_end = max(cur_end, h["norm_end"])
                cur_hits.append(h)
            else:
                loci.append(_make_locus(strain, chrom, strand, cur_start, cur_end, cur_hits))
                cur_start = h["norm_start"]
                cur_end   = h["norm_end"]
                cur_hits  = [h]
        loci.append(_make_locus(strain, chrom, strand, cur_start, cur_end, cur_hits))

    return loci


def _make_locus(strain, chrom, strand, start, end, hits):
    return {
        "strain":     strain,
        "chrom":      chrom,
        "strand":     strand,
        "norm_start": start,
        "norm_end":   end,
        "bitscore":   sum(h["bitscore"] for h in hits),
        "bitscore_max": max(h["bitscore"] for h in hits),
        "n_hits":     len(hits),
        "pident_max": max(h["pident"] for h in hits),
        "evalue_min": min(h["evalue"] for h in hits),
        "length":     max(h["length"] for h in hits),
        "source":     hits[0].get("source", "blast"),
    }


def select_loci(loci, paralog_frac=PARALOG_FRAC):
    """Per strain: keep loci with bitscore_max >= paralog_frac × best."""
    by_strain = defaultdict(list)
    for l in loci:
        by_strain[l["strain"]].append(l)

    selected = {}
    for strain, sloci in by_strain.items():
        best_bs = max(l["bitscore_max"] for l in sloci)
        threshold = best_bs * paralog_frac
        kept = [l for l in sloci if l["bitscore_max"] >= threshold]
        kept.sort(key=lambda x: -x["bitscore_max"])
        selected[strain] = kept

    return selected


def extract_seq(locus, genome_dict):
    chrom    = locus["chrom"]
    if chrom not in genome_dict:
        return None
    chrom_seq = genome_dict[chrom].seq
    chrom_len = len(chrom_seq)
    start = max(0, locus["norm_start"] - 1 - FLANK)
    end   = min(chrom_len, locus["norm_end"] + FLANK)
    subseq = chrom_seq[start:end]
    if locus["strand"] == "-":
        subseq = subseq.reverse_complement()
    return str(subseq)


def load_genome(strain):
    gf = GENOMES_DIR / f"{strain}.fasta"
    return SeqIO.to_dict(SeqIO.parse(gf, "fasta"))


def process_gene(gene):
    blast_dir = BLAST_DIR / gene
    seq_dir   = SEQ_DIR   / gene

    blastn_file  = blast_dir / f"{gene}_blastn.tsv"
    tblastx_file = blast_dir / f"{gene}_tblastx.tsv"

    blastn_hits  = parse_blast(blastn_file)  if blastn_file.exists()  else []
    tblastx_hits = [h for h in (parse_blast(tblastx_file) if tblastx_file.exists() else [])
                    if h["bitscore"] >= MIN_TBLASTX_BS and h["pident"] >= MIN_TBLASTX_PIDENT]

    # Tag source
    for h in blastn_hits:  h["source"] = "blastn"
    for h in tblastx_hits: h["source"] = "tblastx"

    blastn_strains = {h["strain"] for h in blastn_hits}
    # Only use tBLASTx for strains NOT found by BLASTn
    tblastx_new = [h for h in tblastx_hits if h["strain"] not in blastn_strains]

    print(f"\n=== {gene} ===")
    print(f"  BLASTn : {len(blastn_hits)} hits in {len(blastn_strains)} strains")
    print(f"  tBLASTx (new strains only): {len(tblastx_new)} hits")

    # Cluster and select loci
    blastn_loci  = cluster_hits(blastn_hits)
    tblastx_loci = cluster_hits(tblastx_new)
    all_loci     = blastn_loci + tblastx_loci

    selected = select_loci(all_loci)

    all_strains = sorted(p.stem for p in GENOMES_DIR.glob("*.fasta"))
    summary_rows = []
    fasta_records = []
    detail_rows   = []

    for strain in all_strains:
        loci = selected.get(strain, [])
        if not loci:
            summary_rows.append({
                "strain": strain, "present": "absent",
                "n_copies": 0, "best_pident": "NA",
                "best_bitscore_max": "NA", "source": "NA"
            })
            continue

        genome_dict = load_genome(strain)
        written = 0
        for i, locus in enumerate(loci, 1):
            seq = extract_seq(locus, genome_dict)
            if not seq:
                continue
            suffix = f"_copy{i}" if len(loci) > 1 else ""
            seq_id = f"{strain}{suffix}"
            fasta_records.append(f">{seq_id}\n{seq}\n")
            detail_rows.append({
                "strain": strain, "copy": i,
                "chrom": locus["chrom"],
                "start": locus["norm_start"], "end": locus["norm_end"],
                "strand": locus["strand"],
                "n_hits": locus["n_hits"],
                "pident_max": f"{locus['pident_max']:.1f}",
                "bitscore_max": f"{locus['bitscore_max']:.1f}",
                "evalue_min": f"{locus['evalue_min']:.2e}",
                "source": locus["source"],
                "region_length": locus["norm_end"] - locus["norm_start"] + 1,
            })
            written += 1

        best = loci[0]
        summary_rows.append({
            "strain": strain, "present": "present",
            "n_copies": written,
            "best_pident": f"{best['pident_max']:.1f}",
            "best_bitscore_max": f"{best['bitscore_max']:.1f}",
            "source": best["source"],
        })

    # Write outputs
    fasta_out = seq_dir / f"{gene}_hits.fasta"
    with open(fasta_out, "w") as f:
        f.writelines(fasta_records)

    summary_out = seq_dir / f"{gene}_presence_absence.tsv"
    with open(summary_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    detail_out = seq_dir / f"{gene}_blast_hits_detail.tsv"
    if detail_rows:
        with open(detail_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(detail_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(detail_rows)

    present_n = sum(1 for r in summary_rows if r["present"] == "present")
    print(f"  Sequences extracted: {len(fasta_records)}")
    print(f"  Present: {present_n}/{len(all_strains)} strains")
    print(f"  → {fasta_out}")
    print(f"  → {summary_out}")


if __name__ == "__main__":
    for gene in GENES:
        process_gene(gene)
    print("\nDone.")
