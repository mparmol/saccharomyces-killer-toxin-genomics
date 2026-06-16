
from pathlib import Path
import pandas as pd

PROJ   = Path(r"")
DET    = PROJ / "" / ""
META   = PROJ / "" / "" / "genome_metadata.tsv"
DETECT = PROJ / "" / "05_hmm_discovery" / "detection_method_summary.tsv"
OUT    = PROJ / "" / "" / "" / "" / ""
OUT.mkdir(parents=True, exist_ok=True)

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
meta["Species"] = meta["Species"].str.strip()
s2spp = dict(zip(meta["Strain"], meta["Species"]))

# phylogenetic strain order by species (reviewer #138), then strain name
SPP_ORDER = ["Saccharomyces cerevisiae","Saccharomyces paradoxus",
             "Saccharomyces mikatae","Saccharomyces jurei",
             "Saccharomyces kudriavzevii","Saccharomyces arboricola",
             "Saccharomyces uvarum","Saccharomyces chiloensis",
             "Saccharomyces eubayanus"]
spp_rank = {s:i for i,s in enumerate(SPP_ORDER)}

ROMAN = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]
def chr_num(c):
    c = c.replace("chr","")
    return ROMAN.index(c)+1 if c in ROMAN else 99

def load_hits(gene):
    df = pd.read_csv(DET / gene / f"{gene}_blast_hits_detail.tsv", sep="\t")
    df["strain"] = df["strain"].str.strip()
    pa = DET / gene / f"{gene}_presence_absence.tsv"
    if pa.exists():
        dfp = pd.read_csv(pa, sep="\t")
        valid = set(dfp.loc[dfp["present"]=="present","strain"].str.strip())
        df = df[df["strain"].isin(valid)]
    return df

khs1 = load_hits("KHS1")
khr1 = load_hits("KHR1")

# nhmmer novel KHR1 loci
nov = pd.read_csv(DETECT, sep="\t")
nov = nov[(nov["gene"]=="KHR1") & (nov["is_novel"]==True)]
nhmmer_strains = set(nov["strain"].str.strip())

def gene_summary(df, strain):
    g = df[df["strain"]==strain]
    if len(g)==0:
        return "absent", "", 0
    chroms = sorted(g["chrom"].unique(), key=chr_num)
    n_hits = len(g)
    return ";".join(c.replace("chr","") for c in chroms), chroms, n_hits

def khs1_architecture(strain):
    g = khs1[khs1["strain"]==strain]
    if len(g)==0:
        return "absent"
    if len(g)==1:
        return "single locus (1 hit)"
    chroms = list(g["chrom"]); strands=list(g["strand"]); starts=sorted(g["start"])
    if len(set(chroms))>1:
        return f"two chromosomal loci ({'+'.join(c.replace('chr','') for c in sorted(set(chroms), key=chr_num))})"
    # same chromosome
    gap = abs(starts[-1]-starts[0])
    if len(set(strands))>1 and gap < 3000:
        return "single dimeric locus (toxin+antitoxin ORFs, opposite strands)"
    return f"two loci on same chromosome (gap ~{gap/1000:.1f} kb)"

rows = []
strains = sorted(set(meta["Strain"]), key=lambda s:(spp_rank.get(s2spp.get(s,""),99), s))
for s in strains:
    spp = s2spp.get(s, "?")
    kc, kchr, kn = gene_summary(khs1, s)
    rc, rchr, rn = gene_summary(khr1, s)
    arch = khs1_architecture(s)
    rnote = ""
    if s in nhmmer_strains:
        rnote = "second divergent KHR1 locus (nhmmer) on chr IX"
    rows.append({
        "Strain": s,
        "Species": spp.replace("Saccharomyces ","S. "),
        "KHS1_chrom": kc,
        "KHS1_BLAST_hits": kn,
        "KHS1_locus_architecture": arch,
        "KHR1_chrom": rc,
        "KHR1_BLAST_hits": rn,
        "KHR1_notes": rnote,
    })

df = pd.DataFrame(rows)
df.to_csv(OUT / "Supplementary_Table_S1_copy_number.tsv", sep="\t", index=False)
try:
    df.to_excel(OUT / "Supplementary_Table_S1_copy_number.xlsx", index=False)
except Exception as e:
    print("xlsx skip:", e)

# quick console summary of the informative (multi-locus) cases
print("Saved Supplementary Table S1:", len(df), "strains")
print("\n-- KHS1 strains with >1 BLAST hit --")
multi = df[df["KHS1_BLAST_hits"]>1]
for _,r in multi.iterrows():
    print(f"  {r['Strain']:12s} {r['Species']:18s} chr{r['KHS1_chrom']:8s} -> {r['KHS1_locus_architecture']}")
print("\n-- KHR1 strains with second nhmmer locus --")
for _,r in df[df["KHR1_notes"]!=""].iterrows():
    print(f"  {r['Strain']:12s} {r['Species']}")
