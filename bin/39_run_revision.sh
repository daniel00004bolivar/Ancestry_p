#!/usr/bin/env bash
# 39_run_revision.sh — revised autosomal analysis (Phase 2 of the revision memo)
#
# One entry point, run from the project root. Inputs: the per-pool ANGSD
# counts from 05_angsd_analysis.sh, data/hgdp/hgdpGeo.txt.gz and the GRCh37
# reference. Pool sizes: POOL1 = 50, POOL2 = 50, HOSPITAL = 200 people.
# Every random step has a fixed seed (20260925).
#
#   30  oriented SNP universe (HGDP proxies, ancestral/derived coding)   A4 A7
#   31  pool read counts at those sites                                  A5
#   32  common mask + QC summary                                         A6 C1
#   33  count-based ancestry, Delta_HUN, POOL1-POOL2, jackknife          A2 A9
#   34  Pool-seq FST (poolfstat ANOVA), jackknife                        A3
#   35  site-by-site depth matching, then 33 + 34 again                  A10
#   36  sensitivity: transversions, filters, reference panels, blocks   S7 S8
#   37  simulation of the full design                                    A8 C5
#   38  pool-vs-reference FST, common mask (supplement)                  A6
set -euo pipefail
cd "$(dirname "$0")/.."

Rscript bin/30_build_oriented_sites.R
bash    bin/31_extract_pool_counts.sh
Rscript bin/32_build_common_mask.R
Rscript bin/33_ancestry_counts.R --tag main
Rscript bin/34_fst_poolfstat.R  --tag main
Rscript bin/35_depth_matched.R
Rscript bin/33_ancestry_counts.R --table out_global/revision/analysis_table_depthmatched.tsv.gz --tag depthmatched
Rscript bin/34_fst_poolfstat.R  --table out_global/revision/analysis_table_depthmatched.tsv.gz --tag depthmatched
Rscript bin/36_sensitivity.R
Rscript bin/37_simulate_validation.R
Rscript bin/38_fst_references.R
