
import subprocess, re, sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(r"")
BLAST_DIR = PROJ / "results" / "02_blast_detection"
HMM_DIR   = PROJ / "results" / "05_hmm_discovery" / ""
OUT       = PROJ / "results" / "09_tree_comparison" / "cdhit"
OUT.mkdir(parents=True, exist_ok=True)


def translate_longest_orf(nucl_seq):
    """Return the longest protein translated from any of the 6 reading frames."""
    best = ""
    seq = str(nucl_seq).upper().replace("-", "")
    rc = str(Seq(seq).reverse_complement())
    for strand_seq in (seq, rc):
        for frame in range(3):
            s = strand_seq[frame:]
            s = s[:len(s) - len(s)%3]
            prot = str(Seq(s).translate())
            # Take the longest ORF (between two stop codons, or start to stop)
            orfs = re.findall(r'M[^*]*', prot)
            for orf in orfs:
                if len(orf) > len(best):
                    best = orf
            # fallback: longest contiguous non-stop segment
            for seg in prot.split('*'):
                if len(seg) > len(best):
                    best = seg
    return best


def fasta_to_proteins(nucl_fasta, out_fasta, label_suffix=""):
    """Translate nucleotide FASTA to protein FASTA (longest ORF per seq)."""
    records = []
    for rec in SeqIO.parse(nucl_fasta, "fasta"):
        prot = translate_longest_orf(rec.seq)
        if len(prot) < 50:
            continue
        new_id = rec.id + label_suffix
        records.append(f">{new_id}\n{prot}")
    with open(out_fasta, "w") as fh:
        fh.write("\n".join(records) + "\n")
    print(f"  {nucl_fasta.name} -> {len(records)} proteins in {out_fasta.name}")
    return len(records)


def run_cdhit(input_fasta, out_prefix, identity, word_size=None):
    """Run CD-HIT at given identity threshold. Returns cluster count."""
    ws = word_size or (5 if identity >= 0.7 else 4 if identity >= 0.6 else 3)
    cmd = [
        "cd-hit",
        "-i", str(input_fasta),
        "-o", str(out_prefix),
        "-c", str(identity),
        "-n", str(ws),
        "-d", "0",       # unlimited description length in .clstr
        "-T", "8",
        "-M", "4000",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("CD-HIT ERROR:", result.stderr[:500])
        return 0
    clstr_file = Path(str(out_prefix) + ".clstr")
    if clstr_file.exists():
        clusters = len(re.findall(r'^>Cluster', clstr_file.read_text(), re.MULTILINE))
        return clusters
    return 0


def parse_clstr(clstr_file):
    """Parse .clstr file -> dict {seq_id: cluster_id}."""
    mapping = {}
    current_cluster = None
    for line in open(clstr_file):
        line = line.strip()
        if line.startswith(">Cluster"):
            current_cluster = int(line.split()[1])
        elif line and current_cluster is not None:
            m = re.search(r'>([^.]+)\.\.\.', line)
            if m:
                mapping[m.group(1)] = current_cluster
    return mapping


# ── Step 1: Translate nucleotide sequences to protein ────────────────────────
print("\n[1] Translating nucleotide sequences to protein...")

khs1_nucl = BLAST_DIR / "KHS1" / "KHS1_primary_hits.fasta"
khr1_nucl = BLAST_DIR / "KHR1" / "KHR1_clean_hits.fasta"
khr1_hmm  = HMM_DIR / "KHR1_new_candidates_nucl.fasta"

khs1_prot_fa = OUT / "KHS1_proteins.fasta"
khr1_prot_fa = OUT / "KHR1_all_proteins.fasta"

n_khs1 = fasta_to_proteins(khs1_nucl, khs1_prot_fa)
n_khr1_blast = fasta_to_proteins(khr1_nucl, OUT / "KHR1_blast_proteins.fasta")
n_khr1_hmm   = fasta_to_proteins(khr1_hmm,  OUT / "KHR1_hmm_proteins.fasta", label_suffix="")

# Merge KHR1 BLAST + HMM into one file
with open(khr1_prot_fa, "w") as fh:
    for part in (OUT / "KHR1_blast_proteins.fasta", OUT / "KHR1_hmm_proteins.fasta"):
        if part.exists():
            fh.write(part.read_text())
print(f"  KHR1 combined: {n_khr1_blast} blast + {n_khr1_hmm} nhmmer proteins")

# ── Step 2: CD-HIT at multiple thresholds ─────────────────────────────────────
print("\n[2] Running CD-HIT at identity thresholds: 90, 80, 70, 60, 50%...")

THRESHOLDS = [0.90, 0.80, 0.70, 0.60, 0.50]
summary_rows = []

for gene, prot_fa in (("KHS1", khs1_prot_fa), ("KHR1", khr1_prot_fa)):
    for t in THRESHOLDS:
    	pct = int(t * 100)
    	prefix = OUT / f"{gene}_cdhit_{pct}"
    	n_clusters = run_cdhit(prot_fa, prefix, t)
    	print(f"  {gene} @ {pct}% -> {n_clusters} clusters")
    	summary_rows.append({"gene": gene, "identity_pct": pct, "n_clusters": n_clusters,
    	                      "clstr_file": str(prefix) + ".clstr"})

# ── Step 3: Parse clusters at 80% (main working threshold) and report ─────────
print("\n[3] Parsing clusters at 80% identity (main threshold)...")

import csv

with open(OUT / "cluster_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene","identity_pct","n_clusters","clstr_file"], delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

print("\nDone. Cluster counts:")
for row in summary_rows:
    print(f"  {row['gene']} @ {row['identity_pct']}% : {row['n_clusters']} clusters")

print(f"\nOutput in: {OUT}")
