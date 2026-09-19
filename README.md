# Ancestria — genetic ancestry in pool-seq populations from the Colombian Caribbean

Reproducible pipeline for estimating genetic ancestry composition from
pooled-sequencing (pool-seq) data of three cohorts from the Colombian
Caribbean coast (POOL1, POOL2, HOSPITAL), using a pool-seq-aware approach
(ANGSD + NGSadmix + PCAngsd, count-based F<sub>ST</sub>, mitochondrial
haplogroup composition).

**This repository is code only, and contains only the scripts that were
actually used to build the manuscript.** No sequencing data, no results, no
manuscript — see `DATA_AVAILABILITY.md` for what's excluded and why, and
`bin/README_VIGENTE.md` for what else was audited but left out (obsolete
iterations, a complementary method not used in this article, and two
extensions that failed validation).

## Structure

```
bin/                      Pipeline code (see bin/README_VIGENTE.md)
data/panel/                Public 1000 Genomes metadata (reference populations)
env.sh                     Project environment variables
DATA_AVAILABILITY.md       What's excluded and how to obtain/reproduce it
```

Sequencing data, intermediate/final outputs, results, figures, and the
manuscript are all excluded from the repository (`.gitignore`) — they are
generated locally by running the pipeline against the (access-restricted)
raw data.

## Start here: `bin/README_VIGENTE.md`

It lists every script in run order with exactly what table/figure in the
manuscript it produces.

## Pipeline (ANGSD / NGSadmix / PCAngsd), in brief

Full detail, flags, and manuscript cross-references in
`bin/README_VIGENTE.md`. Short version:

1. `05_angsd_analysis.sh` — per-pool allele frequencies, counts, IBS, SFS/F<sub>ST</sub>
2. `05a_fst_from_counts.R`, `05d_fst_per_population.R` — F<sub>ST</sub> tables (inter-pool and vs. reference populations)
3. `06_build_5pop_panel.sh`, `11_full_genome_panels.sh` — reference panel construction from 1000 Genomes
4. `07_build_ref_beagle.R` / `07b_build_ref_beagle_v2.R`, `08_merge_beagle_ngsadmix.R` / `08b_merge_beagle_ngsadmix_v2.R` — reference genotype likelihoods + merge with pools (continental and granular panels)
5. `12_rebuild_and_run.sh` — NGSadmix / PCAngsd / Evanno replicates, full genome
6. `10_evanno_deltaK.R` — optimal K selection
7. `09_plot_admixture_pca.R` / `09b_plot_admixture_pca_v2.R`, `14_pool_focused_plots.R`, `15_final_figures.R` — figures
8. `16_pool_replicates.sh` + `16b_summarize_replicates.R` — uncertainty quantification (block jackknife by chromosome)
9. `17_summarize_mixemt.R` (+ `mixemt` installation) — mitochondrial haplogroup composition
10. `20_new_findings_figures.R` — figures for the jackknife, chromosome X F<sub>ST</sub>, and mtDNA results

## Third-party tools (not vendored — clone at the exact commit used)

`tools/angsd_src/` is excluded from this repository: it is itself a full git
clone of an external project, and vendoring a nested `.git` history inside
this repo causes more problems than it solves. Reconstruct it with:

```bash
# ANGSD (includes NGSadmix and realSFS in misc/), used at commit 6b5d906
git clone https://github.com/ANGSD/angsd.git tools/angsd_src
cd tools/angsd_src && git checkout 6b5d906 && make
```

## Environment and dependencies

Exact versions used and verified in this project:

| Tool | Version | Notes |
|---|---|---|
| ANGSD | 6b5d906 (htslib 1.23.1) | build from [source](https://github.com/ANGSD/angsd); `NGSadmix` and `realSFS` ship in `misc/` |
| PCAngsd | 1.36.4 | `pip install pcangsd` |
| bcftools / samtools | 1.23.1 | |
| R | 4.5.2 | packages: `data.table`, `ggplot2`, `quadprog`, `nnls`, `scales`, `ggrepel`, `patchwork` |
| mixemt | 0.1 | see `tools/mixemt_requirements.txt` — needs a Python venv with `setuptools<81` |
| Python | 3.x | only for the `mixemt` venv |

Base path variables live in `env.sh` (auto-detected relative to the
script's location, portable across machines).

### Installing mixemt

```bash
python3 -m venv tools/mixemt_venv
source tools/mixemt_venv/bin/activate
pip install -r tools/mixemt_requirements.txt
git clone https://github.com/svohr/mixemt.git tools/mixemt
pip install -e tools/mixemt
```
