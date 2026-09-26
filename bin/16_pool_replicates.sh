#!/usr/bin/env bash
# 16_pool_replicates.sh — uncertainty quantification for NGSadmix K=5 via
# block jackknife by chromosome (leave-one-chromosome-out).
#
# CONTEXT / why jackknife by chromosome and not read subsampling:
# subsampling each pool into "individual-depth" replicates (~20-30x) to
# generate independent pseudo-individuals, analogous to iAdmix's 10-run
# bootstrap, was considered first. Discarded: real combined pool depth is
# only ~3-4x (confirmed with `samtools depth` on a test window), i.e.
# ~0.08x per individual within the pool -- no room to subsample without
# running out of data. Instead, this script quantifies uncertainty by
# re-running NGSadmix 22 times, excluding a different chromosome each time,
# on the SAME Beagle already used for the manuscript's results
# (out_global/angsd/ngsadmix/combined_full.beagle.gz, read-only). The
# spread of Q proportions across the 22 runs gives a jackknife standard
# error per pool/component, comparable in spirit to the SD iAdmix reports
# in results/iadmix_summary.tsv.
#
# METHODOLOGICAL WARNING (not resolved by this script): the underlying
# problem remains -- each pool is ONE diploid pseudo-individual in the
# Beagle. The jackknife quantifies how much the estimate varies depending
# on which genome region is used; it does NOT correct the structural bias
# of treating a pool of 50-200 people as 2 chromosome copies. See
# bin/README_VIGENTE.md.
#
# Usage: bash bin/16_pool_replicates.sh [--k 5] [--threads 4]
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ_DIR"

K=5
THREADS=4
while [[ $# -gt 0 ]]; do
    case "$1" in
        --k) K="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

SRC_BEAGLE="out_global/angsd/ngsadmix/combined_full.beagle.gz"
NGSADMIX="$PROJ_DIR/tools/angsd_src/misc/NGSadmix"
OUT_DIR="out_global/angsd/pool_replicates"
WORK_DIR="work/pool_replicates"
mkdir -p "$OUT_DIR" "$WORK_DIR"

if [[ ! -f "$SRC_BEAGLE" ]]; then
    echo "ERROR: $SRC_BEAGLE does not exist (this script only READS the already-built beagle, it does not rebuild it)" >&2
    exit 1
fi

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "=== Chromosome jackknife (leave-one-chromosome-out), K=$K, beagle=$SRC_BEAGLE ==="
log "NOTE: this is an exploratory uncertainty estimate, separate from the results already cited in the manuscript."

# set +o pipefail inside the subshell: `head -1` closes the pipe before
# gzcat finishes writing, which triggers SIGPIPE in gzcat (exit 141) and
# would abort the script under pipefail even though head got the line fine.
HEADER=$(set +o pipefail; gzcat "$SRC_BEAGLE" | head -1)

for CHR in $(seq 1 22); do
    SUBSET="$WORK_DIR/excl_chr${CHR}.beagle"
    log "Excluding chr${CHR}: building subset..."
    {
        echo "$HEADER"
        gzcat "$SRC_BEAGLE" | tail -n +2 | awk -F'\t' -v chr="$CHR" '$1 !~ ("^" chr "_")'
    } > "$SUBSET"
    gzip -f "$SUBSET"

    log "Excluding chr${CHR}: running NGSadmix K=$K..."
    "$NGSADMIX" -likes "${SUBSET}.gz" \
        -K "$K" -P "$THREADS" -seed 1 -minMaf 0.01 \
        -outfiles "$OUT_DIR/jk_excl_chr${CHR}" \
        > "$OUT_DIR/jk_excl_chr${CHR}.runlog" 2>&1

    rm -f "${SUBSET}.gz"
done

log "=== Jackknife complete: 22 runs in $OUT_DIR ==="
log "Next step: Rscript bin/16b_summarize_replicates.R"
