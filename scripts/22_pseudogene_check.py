
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import pandas as pd
import re

PROJ    = Path(r"")
ALN_DIR = PROJ / "re" / ""
META    = PROJ / "re" / "01_metadata" / "genome_metadata.tsv"
OUT     = PROJ / "re" / ""

CORR = {
    "UFRJ50816T":"UFRJ50816","CDFM21L.1":"CDFM21L1","MNFM4L.2":"MNFM4L2",
    "RTSP5L.1":"RTSP5L1","SN-QL4-2":"SNQL42","LL2011-012":"LL2011012",
    "LL2012-001":"LL2012001","LL2012-016":"LL2012016","LL2012-018":"LL2012018",
    "MSH-604":"MSH604","BR6-2":"BR62","JXXY16.1":"JXXY161",
    "UWOPS03-461.4":"UWOPS034614","UWOPS91-917.1":"UWOPS919171","XXYS1.4":"XXYS14",
    "EM14S01-3B":"EM14S013B",
}
meta = pd.read_csv(META, sep="\t")[["Strain","Species"]].copy()
meta["Strain"] = meta["Strain"].str.strip().replace(CORR)
s2spp = dict(zip(meta["Strain"], meta["Species"]))


def check_alignment_pseudogenes(aln_file, gene_name):
    """
    Translate each sequence in the alignment in frame 0 (the CDS frame).
    Flag sequences with in-frame stop codons.

    Note: trimAl alignments keep columns that are present in >=50% of sequences.
    Stop codons in the alignment (after removing gaps) indicate pseudogenisation.
    """
    rows = []
    for rec in SeqIO.parse(aln_file, "fasta"):
        strain = str(rec.id).strip()
        spp = s2spp.get(strain, "?").replace("Saccharomyces ","S. ")

        # Remove gaps for translation
        seq_nogap = str(rec.seq).upper().replace("-", "").replace("N", "A")  # Ns → A as placeholder
        n_len = len(seq_nogap)

        # Try all 3 forward frames, pick the one with fewest stops
        best_frame = 0
        best_prot = ""
        min_stops = 9999
        for frame in range(3):
            s = seq_nogap[frame:]
            s = s[:len(s) - len(s)%3]
            if len(s) < 30:
                continue
            prot = str(Seq(s).translate())
            # Count internal stops (exclude final stop)
            prot_nostop = prot.rstrip("*")
            n_stops = prot_nostop.count("*")
            if n_stops < min_stops:
                min_stops = n_stops
                best_prot = prot
                best_frame = frame

        # ORF analysis
        prot_internal = best_prot.rstrip("*")
        n_internal_stops = prot_internal.count("*")
        n_aa = len(prot_internal)
        has_stop = best_prot.endswith("*")

        flags = []
        if n_internal_stops > 0:
            flags.append(f"internal_stop(n={n_internal_stops})")
        if not has_stop:
            flags.append("no_stop_codon(truncated?)")

        rows.append({
            "strain": strain,
            "species": spp,
            "seq_len_nt": n_len,
            "best_frame": best_frame,
            "n_aa": n_aa,
            "n_internal_stops": n_internal_stops,
            "has_stop_codon": has_stop,
            "flags": "|".join(flags) if flags else "ok",
        })

    return pd.DataFrame(rows)


print("=" * 65)
print("PSEUDOGENE CHECK FROM MAFFT TRIMMED ALIGNMENTS")
print("Note: stop codons expected from frameshifts/divergence in non-CDS")
print("regions — only repeated stops across many sequences are meaningful.")
print("=" * 65)

GENE_FILES = {
    "KHS1": ALN_DIR / "KHS1" / "KHS1_primary_mafft_trim.fasta",
    "KHR1": ALN_DIR / "KHR1" / "KHR1_mafft_gt50.fasta",
}

for gene, aln_file in GENE_FILES.items():
    if not aln_file.exists():
        print(f"\n{gene}: alignment not found at {aln_file}")
        continue

    df = check_alignment_pseudogenes(aln_file, gene)
    out_tsv = OUT / f"pseudogene_check_{gene}.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    n_with_stops = (df["n_internal_stops"] > 0).sum()
    median_aa = df["n_aa"].median()

    print(f"\n── {gene}  (n={len(df)}, alignment file: {aln_file.name}) ──")
    print(f"   Median protein length (in best frame): {median_aa:.0f} aa")
    print(f"   Sequences with ≥1 internal stop codon: {n_with_stops}/{len(df)}")
    print(f"   Sequences without terminal stop (truncated): {(~df['has_stop_codon']).sum()}/{len(df)}")

    # Report only the cleanest sequences (fewest stops = likely functional)
    clean = df[df["n_internal_stops"] == 0].sort_values("n_aa", ascending=False)
    print(f"\n   CLEAN (0 internal stops): {len(clean)} sequences")
    for _, row in clean.iterrows():
        print(f"     {row['strain']:15s} {row['species']:20s} {row['n_aa']:4d} aa  frame={row['best_frame']}")

    # Report highly pseudogenised (≥5 stops)
    many_stops = df[df["n_internal_stops"] >= 5].sort_values("n_internal_stops", ascending=False)
    print(f"\n   HEAVILY PSEUDOGENISED (≥5 internal stops): {len(many_stops)}")
    for _, row in many_stops.head(15).iterrows():
        print(f"     {row['strain']:15s} {row['species']:20s} {row['n_aa']:4d} aa  stops={row['n_internal_stops']}")

    # Report 1-4 stops (possible pseudogenes)
    some_stops = df[(df["n_internal_stops"] >= 1) & (df["n_internal_stops"] < 5)]
    print(f"\n   POSSIBLE PSEUDOGENES (1-4 internal stops): {len(some_stops)}")
    for _, row in some_stops.sort_values("n_internal_stops").iterrows():
        print(f"     {row['strain']:15s} {row['species']:20s} {row['n_aa']:4d} aa  stops={row['n_internal_stops']}")

print(f"\nTSVs saved to: {OUT}")
