#!/usr/bin/env bash
# 31_extract_pool_counts.sh — ACGT read counts of each pool at the oriented sites
#
# Input : out_global/angsd/per_pool/<POOL>_counts.{pos,counts}.gz
#         (ANGSD -doCounts 1 -dumpCounts 3 -minMapQ 30 -minQ 20 -setMinDepth 5;
#          see 05_angsd_analysis.sh)
#         out_global/revision/sites_oriented.tsv.gz (bin/30_build_oriented_sites.R)
# Output: out_global/revision/counts_<POOL>.tsv.gz  (chr pos depth A C G T)
#
# All four base counts are kept, so zero derived-allele counts survive as data
# (memo A5). Sites below ANGSD's depth floor of 5 are simply absent.
set -euo pipefail
cd "$(dirname "$0")/.."

SITES="out_global/revision/sites_oriented.tsv.gz"
IN_DIR="out_global/angsd/per_pool"
OUT_DIR="out_global/revision"

extract() {
    local pool=$1
    paste <(gzcat "$IN_DIR/${pool}_counts.pos.gz") <(gzcat "$IN_DIR/${pool}_counts.counts.gz") |
    awk -F'\t' -v OFS='\t' -v sites=<(gzcat "$SITES") '
        BEGIN { while ((getline line < sites) > 0) { split(line, f, "\t"); if (f[1] != "chr") keep[f[1] ":" f[2]] = 1 }
                print "chr", "pos", "depth", "A", "C", "G", "T" }
        FNR > 1 && (($1 ":" $2) in keep) { print $1, $2, $3, $4, $5, $6, $7 }
    ' | gzip > "$OUT_DIR/counts_${pool}.tsv.gz"
    echo "$pool: $(gzcat "$OUT_DIR/counts_${pool}.tsv.gz" | tail -n +2 | wc -l | tr -d ' ') sites"
}

for pool in POOL1 POOL2 HOSPITAL; do extract "$pool" & done
wait
