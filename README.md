# Ancestry in pooled sequencing libraries from the Colombian Caribbean

Code for estimating continental ancestry from pooled DNA sequencing of three
cohorts from Barranquilla, Colombia — two university-derived libraries of 50
people each (POOL1, POOL2) and one hospital-derived library (HOSPITAL).

A pooled library gives allele counts for a cohort, not genotypes for its
members, so both the estimator and its uncertainty are built around that.
Ancestry is estimated from read counts with a beta-binomial model that takes
the real number of people in each pool; differentiation uses a Pool-seq-aware
F<sub>ST</sub>; and uncertainty combines a block jackknife over SNPs with the
number of people each library effectively represents, a quantity the pipeline
estimates from the data and validates against simulated pools of known
composition.

**Start at [`bin/README_VIGENTE.md`](bin/README_VIGENTE.md)** — it has the flow
diagram, every script with what it produces, and the two silent failure modes
worth knowing about before running anything.

## Running it

```bash
Rscript bin/00_check_requirements.R     # reports anything missing, before a long run starts
bash    bin/01_fetch_reference_data.sh  # downloads the two public reference files
bash    bin/run_pipeline.sh             # the whole analysis
```

Requirements: R ≥ 4.1 with `data.table` (≥ 1.14) and `poolfstat` (≥ 3.0);
samtools; gzip; awk; a C compiler for the two iAdmix cross-checks; and ANGSD,
but only for `bin/05`, which produces the read counts. Every random step is
seeded (20260925).

## What is and is not here

This repository is **code only**. No sequencing data, no derived per-site data,
no results, no manuscript. See [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)
for exactly what can be reproduced with and without access to the alignments,
and for the controlled-access route to the underlying data.

```
bin/                 Pipeline, validation and documentation
data/panel/          Public 1000 Genomes sample metadata
DATA_AVAILABILITY.md What can and cannot be reproduced, and how to request access
```

The reference data the pipeline needs — the UCSC `hgdpGeo` allele-frequency
table and the GRCh37/b37 genome — are public but too large to carry here;
`bin/01_fetch_reference_data.sh` downloads and indexes them.

## Third-party tools

ANGSD is not vendored, since it is itself a git repository. Reconstruct it at
the commit used here:

```bash
git clone https://github.com/ANGSD/angsd.git tools/angsd_src
cd tools/angsd_src && git checkout 6b5d906 && make
```

iAdmix ships as source in `data/iadmix/`; `bin/44` and `bin/45` compile it on
the fly, because the binary distributed with it is Linux x86-64 only.

| Tool | Version used |
|---|---|
| ANGSD | 6b5d906 (htslib 1.23.1) |
| samtools, bcftools | 1.23.1 |
| R | 4.5.2, with `data.table` 1.18.2 and `poolfstat` 3.1.0 |
| bwa | 0.7.17 and 0.7.19 |
| Sentieon | 202010.04 |
| PCAngsd | 1.36.4 (exploratory analyses only) |

## History

The analysis was rebuilt in September 2026 after a methodological review. The
earlier route estimated ancestry with NGSadmix and PCAngsd, which represent a
pool as a single diploid individual and discard the pool size; those results
are retained as exploratory supplementary material and carry none of the
conclusions. Seven scripts whose output was retracted or dropped — the earlier
F<sub>ST</sub> route, the chromosome X and mitochondrial analyses — were removed
from the working tree and remain in the git history. The reasons are set out in
[`bin/README_VIGENTE.md`](bin/README_VIGENTE.md).
