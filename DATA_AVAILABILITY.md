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

Only `bin/` (all pipeline code), `env.sh`, and top-level documentation
(`README.md`, this file, `bin/README_VIGENTE.md`).

## Requesting the raw data

To request access to the underlying BAMs or VCFs for research purposes,
contact the corresponding authors (see the manuscript's *Data Availability*
section — pending completion with the formal controlled-access mechanism
before final publication).

## Reproducing results from scratch

Given access to the BAMs (obtained via the route above) and following
`bin/README_VIGENTE.md` for script order, the full pipeline is reproducible
end to end: the scripts in `bin/` regenerate every intermediate and final
file this repository does not include.
