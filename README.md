# Ancestria — genetic ancestry in pool-seq populations from the Colombian Caribbean

Reproducible pipeline for estimating genetic ancestry composition from
pooled-sequencing (pool-seq) data of three cohorts from the Colombian
Caribbean coast (POOL1, POOL2, HOSPITAL), using a pool-seq-aware approach
(ANGSD + NGSadmix + PCAngsd, count-based F<sub>ST</sub>, mitochondrial
haplogroup composition).

**This repository is code only.** No sequencing data, no results, no
manuscript — see `DATA_AVAILABILITY.md` for what's excluded and why.

## Structure

```
bin/                      All pipeline code (see bin/README_VIGENTE.md)
data/panel/                Public 1000 Genomes metadata (reference populations)
data/iadmix/                Not vendored — see "Third-party tools" below
env.sh                     Project environment variables
DATA_AVAILABILITY.md       What's excluded and how to obtain/reproduce it
```

Sequencing data, intermediate/final outputs, results, figures, and the
manuscript are all excluded from the repository (`.gitignore`) — they are
generated locally by running the pipeline against the (access-restricted)
raw data.

## Start here: `bin/README_VIGENTE.md`

`bin/` accumulated ~40 scripts over the course of development. **Always
start with `bin/README_VIGENTE.md`** — it lists exactly which script
corresponds to each pipeline step, which ones were archived to
`bin/archive/` as obsolete (do not use), and which exploratory results
failed validation (do not cite).

## Main pipeline (ANGSD / NGSadmix / PCAngsd)

General order (full detail, flags, and outputs in `bin/README_VIGENTE.md`):

1. `05_angsd_analysis.sh` — per-pool allele frequencies, counts, IBS, SFS/F<sub>ST</sub> (autosomes or `--chr X`)
2. `05a_fst_from_counts.R` — Hudson-Bhatia F<sub>ST</sub> from nucleotide counts (the robust estimator for pool-seq)
3. `06_build_5pop_panel.sh`, `11_full_genome_panels.sh` — reference panel construction from 1000 Genomes
4. `07b_build_ref_beagle_v2.R`, `08b_merge_beagle_ngsadmix_v2.R` — reference genotype likelihoods + merge with pools
5. NGSadmix / PCAngsd — admixture and principal components
6. `10_evanno_deltaK.R` — optimal K selection (Evanno method)
7. `16_pool_replicates.sh` + `16b_summarize_replicates.R` — uncertainty quantification (block jackknife by chromosome)
8. `17_summarize_mixemt.R` (+ `mixemt` installation) — mitochondrial haplogroup composition
9. `18_build_X_panel_and_run.sh`, `19a_prepare_y_markers_and_counts.sh` + `19_y_haplogroup_mixture.R` — exploratory extensions to chromosomes X and Y (**not validated** — see the warning in each script's header)

## Third-party tools (not vendored — clone at the exact commit used)

`tools/angsd_src/` and `data/iadmix/` are excluded from this repository:
both are themselves full git clones of external projects, and vendoring a
nested `.git` history inside this repo causes more problems than it solves.
Reconstruct them with:

```bash
# ANGSD (includes NGSadmix and realSFS in misc/), used at commit 6b5d906
git clone https://github.com/ANGSD/angsd.git tools/angsd_src
cd tools/angsd_src && git checkout 6b5d906 && make

# iAdmix (complementary pool-seq-aware ancestry method), used at commit c188a70
git clone https://github.com/eliorav/iAdmix.git data/iadmix
cd data/iadmix && git checkout c188a70 && make
```

`data/iadmix/` was not audited in depth as part of this review — see the
"Out of scope" section in `bin/README_VIGENTE.md`.

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

### Installing mixemt (step 8)

```bash
python3 -m venv tools/mixemt_venv
source tools/mixemt_venv/bin/activate
pip install -r tools/mixemt_requirements.txt
git clone https://github.com/svohr/mixemt.git tools/mixemt
pip install -e tools/mixemt
```

## Reproducibility

Every script written or corrected during this audit documents in its header
what it does, what it depends on, and — where relevant — why an exploratory
result should not be treated as valid (see `bin/18_build_X_panel_and_run.sh`
and `bin/19a_prepare_y_markers_and_counts.sh` in particular).
