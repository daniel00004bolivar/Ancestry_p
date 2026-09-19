# Current scripts, by pipeline step

This file exists because `bin/` accumulated ~30 obsolete/duplicate
iterations (`_v2`, `_FIXED`, `_BACKUP`, `_ORIGINAL`, `(copia)`) with no clear
indication of which version was actually in use. Those versions were moved
to `archive/` (not deleted). This table says which one to run at each step.

## ANGSD / NGSadmix / PCAngsd pipeline (audited, current)

| Step | Current script | What it does |
|---|---|---|
| 05 | `05_angsd_analysis.sh` | Genotype likelihoods, ACGT counts, IBS, SAF/SFS per pool |
| 05a | `05a_fst_from_counts.R` | Hudson-Bhatia F<sub>ST</sub> from counts (pool-safe, no diploid assumption) |
| 05b | `05b_compare_angsd_iadmix.R` | Compares ANGSD vs. iAdmix results |
| 05c | `05c_continental_ancestry_v3.R` | Continental AIMs from counts (the only version cited in WORKFLOW.md) |
| 05d | `05d_fst_per_population.R` | F<sub>ST</sub> per reference population |
| 06 | `06_build_5pop_panel.sh` | Builds a 5-superpopulation reference panel from 1000G |
| 07 | `07b_build_ref_beagle_v2.R` | Builds the reference Beagle (v2 = current, supersedes `07_build_ref_beagle.R`) |
| 08 | `08b_merge_beagle_ngsadmix_v2.R` | Merges pools + reference Beagle, runs NGSadmix (v2 = current) |
| 09 | `09b_plot_admixture_pca_v2.R` | PCAngsd + admixture plots (v2 = current) |
| 10 | `10_evanno_deltaK.R` | Evanno method for choosing K |
| 11 | `11_full_genome_panels.sh` | Full-genome panels |
| 12 | `12_rebuild_and_run.sh` | Full rebuild (portability fixed in this audit) |
| 13 | `13_pca_3d.R` | 3D PCA |
| 14 | `14_pool_focused_plots.R` | Pool-focused plots |
| 15 | `15_final_figures.R` | Final manuscript figures |
| 16 | `16_pool_replicates.sh` + `16b_summarize_replicates.R` | **New** — uncertainty quantification via chromosome block jackknife (see `results/pool_jackknife_summary.tsv`) |
| — | `05_angsd_analysis.sh --chr X` | **New, validated** — F<sub>ST</sub>/counts on chromosome X (see `out_global/angsd_X/fst/`) |
| 18 | `18_build_X_panel_and_run.sh` | **New, NOT validated** — full NGSadmix/PCAngsd pipeline on X. See status detail below. Do not use these numbers. |
| 17 | `17_summarize_mixemt.R` | **New, complete** — mtDNA haplogroups via `mixemt` (see `results/mtdna_haplogroups.tsv`) |
| 19a | `19a_prepare_y_markers_and_counts.sh` | **New, NOT validated** — downloads Y VCF, selects markers, ANGSD counts per pool |
| 19 | `19_y_haplogroup_mixture.R` | **New, NOT validated** — Y-chromosome mixture deconvolution (needs 19a first). See `docs/Y_HAPLOGROUP_STATUS.md`. Do not use these numbers. |

### Detail: NGSadmix/PCAngsd on chromosome X (exploratory, not validated)

A pipeline extending NGSadmix/PCAngsd to chromosome X was built from
scratch, since the autosomal HapMap3 panel has no X coverage. Reproducible
with **`bin/18_build_X_panel_and_run.sh`** (requires
`bin/05_angsd_analysis.sh --chr X` to have run first). Produces:
- `work/angsd_X_panel/sites_X.txt` — de novo sites from ANGSD's per-pool MAF calls
- `out_global/angsd_X/ngsadmix/pools_X.beagle.gz` + `refX.beagle` + `combined_X.beagle.gz`
- `out_global/angsd_X/ngsadmix/results/K5_X.qopt`, `out_global/angsd_X/pcangsd/pca_X.cov`

The numeric result **fails a basic sanity check** (all three pools come out
nearly identical, contradicting chromosome X's own F<sub>ST</sub> and every
other line of evidence) and PCAngsd did not converge. Full detail and likely
cause in `docs/NGSADMIX_X_STATUS.md`. The pipeline is reusable
for a retry with better marker selection, but the current numbers **are
not** in any results document or in the manuscript.

## ⚠️ Methodological limitation: pools treated as diploid pseudo-individuals

`NGSadmix` and `PCAngsd` consume the Beagle format, which encodes only 3
states per site (AA/Aa/aa) — i.e., it is diploid by design. Each pool
(POOL1/POOL2/HOSPITAL, ~50 real individuals each) is represented as **a
single diploid pseudo-individual**, not as a weighted population frequency.
This is not a bug fixable with a flag: there is no way to tell
NGSadmix/PCAngsd "this is a pool of 50 people."

- Analyses based on direct ACGT counts (`05a_fst_from_counts.R`,
  `05c_continental_ancestry_v3.R`) **do not have this problem** — they
  compute allele frequency without ever going through a genotype.
- `05_angsd_analysis.sh --dumpCounts 3 -doMajorMinor 4` (the counts block)
  also avoids the problem.
- `bin/16_pool_replicates.sh` quantifies uncertainty via a delete-one-
  chromosome block jackknife (22 replicates) on the same Beagle used for the
  published K=5 result. It does **not** correct the underlying pseudo-
  individual bias — it only measures how much the estimate varies with
  which part of the genome is used. Read-subsampling to build individual-
  depth pseudo-replicates (analogous to iAdmix's 10-run bootstrap) was
  attempted first and abandoned: real combined pool depth is only ~3-4x, far
  too low to subsample into individual-depth replicates.

See the header comment in each script listed above for the full warning
text.

## bin/archive/

Contains the ~30 obsolete/duplicate versions moved during this cleanup
(panel-building scripts with different pool sizes/filters, backups,
copies). Kept intact for traceability, but should not be used to generate
new results — use the table above.

## Out of scope

The `iAdmix` scripts (`iadmix_master_script.sh`, `run_iadmix_*.sh`,
`03_iadmix_genome.sh`, etc.) were not audited in depth in this work — they
were explicitly excluded from scope. Only their obvious duplicates were
archived as general housekeeping.
