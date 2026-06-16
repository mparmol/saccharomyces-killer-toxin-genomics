# KHS1 and KHR1 killer toxin genes in 83 *Saccharomyces* genomes — Analysis pipeline

**Pármol M, Peris F, Peris D** (manuscript in preparation)

This repository contains the analysis scripts used to detect, characterise, and phylogenetically analyse two killer-toxin loci, **KHS1** (NUPAV/YSC0044) and **KHR1** (YSC0002), across 83 assembled nuclear genomes of the genus *Saccharomyces*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
*(DOI will be populated after Zenodo release)*

---

## Software requirements

| Tool | Version | Install |
|------|---------|---------|
| Python | ≥ 3.10 | conda |
| Biopython | ≥ 1.85 | `conda install -c conda-forge biopython` |
| pandas | ≥ 1.5 | `conda install pandas` |
| matplotlib / seaborn | any recent | `conda install matplotlib seaborn` |
| BLAST+ | ≥ 2.16 | `conda install -c bioconda blast` |
| MAFFT | ≥ 7.520 | `conda install -c bioconda mafft` |
| trimAl | ≥ 1.4 | `conda install -c bioconda trimal` |
| IQ-TREE 2 | ≥ 2.3 | `conda install -c bioconda iqtree` |
| HMMER | ≥ 3.4 | `conda install -c bioconda hmmer` |
| CD-HIT | ≥ 4.8.1 | `conda install -c bioconda cd-hit` |
| BUSCO | ≥ 5.7 | `conda install -c bioconda busco` |
| Entrez Direct | any | `conda install -c bioconda entrez-direct` |
| DendroPy | ≥ 4.6 | `conda install -c conda-forge dendropy` |
| python-docx | any | `pip install python-docx` |

A ready-to-use conda environment can be created with:

```bash
conda create -n sac_killer python=3.10 -y
conda activate sac_killer
conda install -c bioconda -c conda-forge \
    biopython pandas matplotlib seaborn scipy openpyxl \
    blast mafft trimal iqtree hmmer cd-hit busco \
    entrez-direct dendropy -y
pip install python-docx
```

---

## Input data

The 83 assembled nuclear genomes are **not included** in this repository (file sizes prohibitive). They are available from the companion study:

> Peris D et al. 2023. Macroevolutionary diversity of traits and genomes in the model yeast genus *Saccharomyces*. *Nat Commun* **14**:690. https://doi.org/10.1038/s41467-023-36354-x

Genome files should be placed in `Nuclear_genomes/` (one FASTA per strain, named `{STRAIN}.fasta`). Headers must follow the format `{STRAIN}__{chr}` (e.g. `CBS432__chrIV`).

The KHS1 (YSC0044) and KHR1 (YSC0002) query sequences used as BLAST/HMM seeds can be downloaded from SGD (https://www.yeastgenome.org/).

---

## Pipeline

Scripts should be run in order. Scripts 01–10 run on Windows (Python) or WSL2 (bash). Scripts 11+ require WSL2 with the `sac_killer` conda environment.

| # | Script | What it does | Engine |
|---|--------|-------------|--------|
| 01 | `01_prepare_blast_db.py` | Merge 83 genomes → single BLAST database | Python/WSL |
| 02 | `02_run_blast.py` | BLASTn + tBLASTx for KHS1 and KHR1 | Python/WSL |
| 03 | `03_extract_sequences.py` | Extract hits, cluster, presence/absence tables | Python/WSL |
| 04 | `04_align_and_tree.sh` | MAFFT + trimAl + IQ-TREE (initial trees) | bash/WSL |
| 05 | `05_figures.py` | Heatmaps, identity distributions, alignment stats | Python |
| 06 | `06_khs1_primary_tree.sh` | IQ-TREE on KHS1 primary copies only (74 seqs) | bash/WSL |
| 07 | `07_tree_figures.py` | Publication-quality tree figures | Python |
| 08 | `08_extra_figures.py` | Chromosomal location, ML heatmap, barplot | Python |
| 09 | `09_khr1_groups_v2.py` | KHR1 group classification (C1/C2/C3) + annotated figure | Python |
| 10 | `10_fix_khr1_full.sh` | Master script: reprocess KHR1 from scratch (false-positive removal) | bash/WSL |
| 12 | `12_copy_number.py` | Copy-number heatmap per strain × chromosome | Python |
| 13 | `13_busco_phylogenomics.sh` + `13b_build_supermatrix.py` | BUSCO whole-genome phylogenomics | bash/WSL |
| 14 | `14_hmm_discovery.sh` + `14b_process_hits.py` | nhmmer HMM discovery of novel homologs | bash/WSL |
| 15 | `15_compare_trees.py` | Gene-tree vs species-tree comparison (RF, tanglegram) | Python/WSL |
| 16b | `16b_nj_boxes.py` | NJ network with bootstrap boxes | Python |
| 16c | `16c_colored_leaves.py` | Tree with colored leaves by species | Python |
| 16e | `16e_splitspy_network.py` | SplitsTree-style network | Python |
| 17 | `17_networks_v2.py` | NeighborNet networks v2 | Python |
| 18 | `18_phylogenomics_tree.py` | Phylogenomics tree figure | Python |
| 19 | `19_copy_number_by_species.py` | Copy number aggregated by species | Python |
| 19b | `19b_copy_heatmap_corrected.py` | Corrected copy heatmap | Python |
| 20 | `20_supp_table_copies.py` | Supplementary table of copy numbers | Python |
| 21 | `21_cdhit_clustering.py` | Protein clustering at 5 identity thresholds (CD-HIT) | Python/WSL |
| 21b | `21b_cdhit_parse_clusters.py` | Parse CD-HIT clusters → strain/chrom annotation | Python/WSL |
| 22 | `22_pseudogene_check.py` | In-frame stop codon detection (pseudogene screen) | Python/WSL |
| 23 | `23_mkiller_hmm.sh` | M-killer dsRNA preprotoxin HMM search in 83 genomes | bash/WSL |

---

## Output structure

```
resultados/
├── 01_metadata/            genome metadata TSV
├── 02_deteccion_blast/     BLAST hits, presence/absence, sequences
├── 03_alineamientos/       MAFFT + trimAl alignments
├── 04_arboles/             IQ-TREE outputs (gene trees)
├── 05_figuras/             All publication figures (PDF + PNG)
├── 06_busco/               BUSCO results per strain
├── 07_phylogenomics/       Supermatrix + species tree
├── 08_discovery/           nhmmer HMM hits + novel candidates
├── 09_comparacion_arboles/ CD-HIT clusters, BLASTp annotations, pseudogene check
└── 10_mkiller_hmm/         M-killer HMM profile + search results
```

---

## Key results summary

| Gene | Present / 83 strains | Key finding |
|------|----------------------|------------|
| KHS1 | 74 (89.2%) | Subtelomeric, highly mobile; found on 11 chromosomes; dimeric NUPAV locus (toxin + antitoxin ORFs in opposite orientation) |
| KHR1 | 68 (81.9%) | Primarily chrIX (C1, 50 strains); chrIII translocations (C2, 16 strains); 2 singletons (C3); 7 divergent loci in *S. kudriavzevii*/*S. jurei* (nhmmer) |

Gene-tree / species-tree discordance (normalised RF: KHS1 = 0.63, KHR1 = 0.70) is consistent with viral capture, horizontal transfer, and subtelomeric recombination rather than vertical descent.

---

## Citation

If you use these scripts, please cite:

> Pármol M, Peris F, Peris D. KHS1 and KHR1 killer-toxin loci in 83 *Saccharomyces* genomes: distribution, chromosomal mobility, and phylogenetic discordance. *[journal, year, DOI — TBD]*

---

## License

Scripts are released under the [MIT License](LICENSE).

The genome data belong to the respective depositors; please follow the data-use terms of the companion study (Peris et al. 2023, Nat Commun).
