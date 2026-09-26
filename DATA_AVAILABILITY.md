# Data availability

This repository is **code only**. Neither raw/derived sequencing data nor
results, figures, or the manuscript are included (see `.gitignore`) — they
are generated locally by running the pipeline.

## Excluded, and why

| What | Would live in | Why it's not here |
|---|---|---|
| POOL1, POOL2, HOSPITAL BAMs | `bam/` | Real participant/patient sequencing data — identifiable, requires a controlled-access repository (dbGaP/EGA) and ethics approval, not a public repo. |
| Pool-level variant VCFs | `entregable_VCF/vcf/`, `out_global/varscan_genome/` | Genomic data derived from the pools at variant level — same reason as above. |
| Genotype likelihoods / per-site counts (Beagle, .mafs, .counts) | `out_global/` | Intermediate data derived from the real samples at individual-site resolution. |
| GRCh37/b37 reference genome | `data/ref/` | Public but heavy (~3GB); downloaded separately (see `bin/README_PIPELINE.md`). |
| 1000 Genomes VCF panels | fetched remotely (streaming) | Public; queried via remote tabix/bcftools streaming, never stored locally in full. |
| Results, figures, tables | `results/`, `manuscrito_G3/`, `entregable_VCF/` | Project decision: this repository ships code only. |

## What IS in this repository

Only `bin/` (the scripts that built the manuscript — see
`bin/README_VIGENTE.md`), `env.sh`, `data/panel/` (public 1000 Genomes
sample metadata), and top-level documentation (`README.md`, this file).

## Requesting the raw data

To request access to the underlying BAMs or VCFs for research purposes,
contact the corresponding authors (see the manuscript's *Data Availability*
section — pending completion with the formal controlled-access mechanism
before final publication).

## What can and cannot be reproduced

Stated precisely, rather than as a blanket claim of full reproducibility.

**Reproducible from this repository alone**, with no access to participant
data: the synthetic test set written by `bin/37_simulate_validation.R`
(`out_global/revision/synthetic_test/`) exercises allele harmonization,
filtering, the direction of the estimator and the summary calculations
against known expected values.

**Reproducible given access to the BAMs** (see the route above), by running
`bash bin/39_run_revision.sh`: the whole revised autosomal analysis —
oriented SNP set, common mask, count-based ancestry model, Pool-seq
F<sub>ST</sub>, depth matching, sensitivity analyses and simulation. Every
random step is seeded (20260925). Reference data are downloaded separately:
the GRCh37 reference genome and the UCSC `hgdpGeo` table of HGDP allele
frequencies (`data/hgdp/`).

**Not reproducible as-is:** the earlier NGSadmix/PCAngsd results
(`bin/07`–`bin/15`), now reported as exploratory. Two of their inputs were
written to `/tmp` during the original run and no longer exist: the site list
`angsd_sites_hapmap3.txt` and the pooled Beagle file
`pools_hapmap3.beagle.gz`. The exact ANGSD command that produced the latter
is preserved in its `.arg` file, but no script in `bin/` regenerates it.

**Not included in any form:** the capture-kit target intervals (BED), which
were never available to this project — see the manuscript's Methods.
