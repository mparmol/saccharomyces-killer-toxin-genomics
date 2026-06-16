#!/usr/bin/env python3


import argparse
import re
import subprocess
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--project_dir", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def parse_nhmmer_tbl(tbl_path):

    hits = []
    for line in Path(tbl_path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 15:
            continue
        hits.append({
            "target":   parts[0],  
            "query":    parts[2],   
            "evalue":   float(parts[12]),
            "score":    float(parts[13]),
            "bias":     float(parts[14]),
            "hmm_from": int(parts[4]),
            "hmm_to":   int(parts[5]),
            "ali_from": int(parts[6]),
            "ali_to":   int(parts[7]),
            "strand":   parts[11],  # '+' o '-'
        })
    return hits


def extract_region(fasta_path, chrom, start, end, strand, flank=0):

    seq = {}
    header = None
    buf = []
    for line in Path(fasta_path).read_text().splitlines():
        if line.startswith(">"):
            if header:
                seq[header] = "".join(buf)
            header = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip().upper())
    if header:
        seq[header] = "".join(buf)

    if chrom not in seq:
        return None

    chrom_seq = seq[chrom]
    s = max(0, start - 1 - flank)
    e = min(len(chrom_seq), end + flank)
    region = chrom_seq[s:e]

    if strand == "-":
        comp = str.maketrans("ACGT", "TGCA")
        region = region.translate(comp)[::-1]
    return region


CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate(seq):
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

def best_orf(nucl_seq, min_len=60):

    comp = str.maketrans("ACGT", "TGCA")
    rev = nucl_seq.translate(comp)[::-1]
    frames = [nucl_seq[i:] for i in range(3)] + [rev[i:] for i in range(3)]
    best = ""
    for frame in frames:
        prot = translate(frame)
        parts = prot.split("*")
        for part in parts:
            if len(part) >= min_len and len(part) > len(best):
                best = part
    return best if best else None


def load_known_coords(project_dir):

    known = {}
    for gene in ["KHS1", "KHR1"]:
        detail = project_dir / "results" / "sequences" / gene / f"{gene}_blast_hits_detail.tsv"
        if not detail.exists():
            continue
        df = pd.read_csv(detail, sep="\t")
        for _, row in df.iterrows():
            key = (str(row["strain"]), str(row["chrom"]))
            known.setdefault(key, []).append((int(row["start"]), int(row["end"])))
    return known


def overlaps(s1, e1, s2, e2, min_overlap=0.5):
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    return overlap / len1 >= min_overlap if len1 > 0 else False


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    project = args.project_dir
    disc_dir = project / "results_2" / "discovery"
    genome_dir = project / "Nuclear genomes"

    known_coords = load_known_coords(project)

    meta = pd.read_csv(project / "results" / "genome_metadata.tsv", sep="\t")
    strain_names = {row["Strain"]: row["Strain"] for _, row in meta.iterrows()}

    ref_lengths = {"KHS1": 2127, "KHR1": 891}

    detect_rows = []

    for gene in ["KHS1", "KHR1"]:
        detail = project / "results" / "sequences" / gene / f"{gene}_blast_hits_detail.tsv"
        if not detail.exists():
            continue
        df = pd.read_csv(detail, sep="\t")
        for _, row in df.iterrows():
            detect_rows.append({
                "strain":           str(row["strain"]),
                "gene":             gene,
                "chrom":            str(row["chrom"]),
                "start":            int(row.get("start", 0)),
                "end":              int(row.get("end", 0)),
                "detection_method": "BLASTn" if gene == "KHS1" else
                                    ("tBLASTx" if row.get("source", "") == "tblastx" else "BLASTn"),
                "pident":           float(row.get("pident_max", 0)),
                "evalue":           float(row.get("evalue_min", 1)),
                "hmm_score":        None,
                "pfam_domain":      None,
                "eggnog_function":  None,
                "is_novel":         False,
                "notes":            "BLAST",
            })

    for gene in ["KHS1", "KHR1"]:
        hits_dir = disc_dir / "hits" / gene
        if not hits_dir.exists():
            continue

        ref_len = ref_lengths.get(gene, 1000)
        new_nucl_seqs = []
        new_prot_seqs = []

        print(f"\n[{gene}] ...")
        n_files = 0
        n_new = 0

        for tbl in sorted(hits_dir.glob("*.tbl")):
            strain = tbl.stem
            n_files += 1
            hits = parse_nhmmer_tbl(tbl)

            fasta = genome_dir / f"{strain}.fasta"
            if not fasta.exists():
                continue

            for hit in hits:
                chrom = hit["target"]
                s, e = hit["ali_from"], hit["ali_to"]
                if s > e:
                    s, e = e, s

                key = (strain, chrom)
                already_known = False
                for ks, ke in known_coords.get(key, []):
                    if overlaps(s, e, ks, ke, min_overlap=0.3):
                        already_known = True
                        break
                if already_known:
                    continue

                region_len = e - s
                if region_len < ref_len * 0.4:
                    continue  
                if hit["evalue"] > 0.001:
                    continue  

                nucl = extract_region(fasta, chrom, s, e, hit["strand"], flank=150)
                if nucl is None or len(nucl) < 100:
                    continue

                prot = best_orf(nucl, min_len=50)

                seq_id = f"{strain}__{chrom}_{s}_{e}_nhmmer"
                new_nucl_seqs.append((seq_id, nucl))
                if prot:
                    new_prot_seqs.append((seq_id, prot))

                detect_rows.append({
                    "strain":           strain,
                    "gene":             gene,
                    "chrom":            chrom,
                    "start":            s,
                    "end":              e,
                    "detection_method": "nhmmer",
                    "pident":           None,
                    "evalue":           hit["evalue"],
                    "hmm_score":        hit["score"],
                    "pfam_domain":      None,
                    "eggnog_function":  None,
                    "is_novel":         True,
                    "notes":            f"nhmmer score={hit['score']:.1f}",
                })
                n_new += 1

        out_nucl = disc_dir / "hits" / f"{gene}_new_candidates_nucl.fasta"
        out_prot = disc_dir / "hits" / f"{gene}_new_candidates_protein.fasta"
        with open(out_nucl, "w") as f:
            for sid, seq in new_nucl_seqs:
                f.write(f">{sid}\n{seq}\n")
        with open(out_prot, "w") as f:
            for sid, seq in new_prot_seqs:
                f.write(f">{sid}\n{seq}\n")

        known_prot = project / "results" / "sequences" / gene
        out_all = disc_dir / "hits" / f"{gene}_all_proteins.fasta"
        with open(out_all, "w") as fout:
            known_fasta = known_prot / (
                "KHS1_primary_hits.fasta" if gene == "KHS1" else "KHR1_clean_hits.fasta"
            )
            if known_fasta.exists():
                for line in known_fasta.read_text().splitlines():
                    if line.startswith(">"):
                        header = line[1:].split()[0]
                        current_nucl = []
                        fout.write(f">{header}_known\n")
                    else:
                        current_nucl.append(line.strip().upper())
                        prot = best_orf("".join(current_nucl), min_len=30)
                        if prot and "\n" not in line:  
                            pass  
            for sid, seq in new_prot_seqs:
                fout.write(f">{sid}\n{seq}\n")

        if n_new == 0:
            print(f"   {gene}")
            print(f"     {gene}")

    df_detect = pd.DataFrame(detect_rows)
    out_detect = disc_dir / "detection_method_summary.tsv"
    df_detect.to_csv(out_detect, sep="\t", index=False)
    print(f"\n: {out_detect}")
    print(f"  : {len(df_detect)}")
    print(df_detect.groupby(["gene", "detection_method", "is_novel"])
                   .size().rename("count").to_string())


if __name__ == "__main__":
    main()
