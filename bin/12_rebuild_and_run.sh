#!/usr/bin/env bash
set -euo pipefail

# METHODOLOGICAL WARNING: NGSadmix and PCAngsd (Step 3 and 4 of this script)
# receive each pool (POOL1/POOL2/HOSPITAL, ~50 real individuals each) as ONE
# diploid pseudo-individual. The Beagle format only allows 3 states per site
# (AA/Aa/aa) -- a structural limitation of the format, not configurable. The
# resulting proportions should be interpreted as the pool's population
# average, not individual ancestry.
# See bin/16_pool_replicates.sh for an uncertainty estimate via block
# jackknife, and bin/README_VIGENTE.md for full detail.

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ_DIR"

WORK_DIR="$PROJ_DIR/work/rebuild_and_run"
mkdir -p "$WORK_DIR"

ANGSD="$PROJ_DIR/tools/angsd_src/angsd"
NGSADMIX="$PROJ_DIR/tools/angsd_src/misc/NGSadmix"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Step 1: Rebuild ref100.beagle and ref110.beagle with all 22 chromosomes ──
log "=== Rebuilding reference Beagle files with all 22 chromosomes ==="

Rscript "$PROJ_DIR/bin/lib_build_ref_beagle_panel.R" \
    "$WORK_DIR/angsd_sites_hapmap3.txt" \
    "$WORK_DIR/ref_extract" \
    "$WORK_DIR/ref_sample_order.txt" \
    "out_global/angsd/ngsadmix/ref100_full.beagle" \
    "v1 (100 individuals)"
log "Panel v1 (100 individuals) rebuilt"

Rscript "$PROJ_DIR/bin/lib_build_ref_beagle_panel.R" \
    "$WORK_DIR/angsd_sites_hapmap3.txt" \
    "$WORK_DIR/ref_extract_v2" \
    "$WORK_DIR/ref_sample_order_v2.txt" \
    "out_global/angsd/ngsadmix/ref110_full.beagle" \
    "v2 (110 individuals)"
log "Panel v2 (110 individuals) rebuilt"

# ── Step 2: Merge with pools beagle ──────────────────────────────────────────
log "=== Merging panels with pools beagle ==="

Rscript -e '
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
WORK_DIR <- args[1]
POOL_NAMES <- c("POOL1","POOL2","HOSPITAL")

# v1: 100 individuals
pools <- fread("out_global/angsd/ngsadmix/pools_hapmap3.beagle.gz", header=TRUE)
ind_idx <- rep(seq_along(POOL_NAMES), each=3)
new_names <- paste0(POOL_NAMES[ind_idx], c("_p1","_p2","_p3"))
setnames(pools, 4:ncol(pools), new_names)

ref <- fread("out_global/angsd/ngsadmix/ref100_full.beagle", header=FALSE)
setnames(ref, 1:3, c("marker","allele1","allele2"))

cat(sprintf("v1: Pools=%d sites, Ref=%d sites\n", nrow(pools), nrow(ref)))
m <- merge(pools, ref, by=c("marker","allele1","allele2"))
cat(sprintf("v1: Combined=%d sites, %d individuals\n", nrow(m), (ncol(m)-3)/3))
fwrite(m, "out_global/angsd/ngsadmix/combined_full.beagle.gz", sep="\t")

# Individuals table
samples <- readLines(file.path(WORK_DIR, "ref_sample_order.txt"))
panel <- fread("data/panel/integrated_call_samples_v3.20130502.ALL.panel")
pop_map <- panel[match(samples, sample), .(sample, pop, super_pop)]
ind_table <- rbind(
    data.table(sample=POOL_NAMES, pop=POOL_NAMES, super_pop="POOL"),
    pop_map
)
fwrite(ind_table, "out_global/angsd/ngsadmix/individuals_full.tsv", sep="\t")

# v2: 110 individuals
ref2 <- fread("out_global/angsd/ngsadmix/ref110_full.beagle", header=FALSE)
setnames(ref2, 1:3, c("marker","allele1","allele2"))

cat(sprintf("v2: Pools=%d sites, Ref=%d sites\n", nrow(pools), nrow(ref2)))
m2 <- merge(pools, ref2, by=c("marker","allele1","allele2"))
cat(sprintf("v2: Combined=%d sites, %d individuals\n", nrow(m2), (ncol(m2)-3)/3))
fwrite(m2, "out_global/angsd/ngsadmix/combined_v2_full.beagle.gz", sep="\t")

samples2 <- readLines(file.path(WORK_DIR, "ref_sample_order_v2.txt"))
pop_map2 <- panel[match(samples2, sample), .(sample, pop, super_pop)]
ind_table2 <- rbind(
    data.table(sample=POOL_NAMES, pop=POOL_NAMES, super_pop="POOL"),
    pop_map2
)
fwrite(ind_table2, "out_global/angsd/ngsadmix/individuals_v2_full.tsv", sep="\t")
cat("Merge complete.\n")
' "$WORK_DIR"

log "Merge complete"

# ── Step 3: Run NGSadmix ─────────────────────────────────────────────────────
log "=== Running NGSadmix (single runs for plotting) ==="
log "NOTE: each pool = 1 diploid pseudo-individual in the Beagle (see warning above)"

RESULTS_V1="out_global/angsd/ngsadmix/results_full"
RESULTS_V2="out_global/angsd/ngsadmix/results_v2_full"
mkdir -p "$RESULTS_V1" "$RESULTS_V2"

# v1: K=2..5
for K in 2 3 4 5; do
    log "NGSadmix v1: K=$K"
    "$NGSADMIX" -likes out_global/angsd/ngsadmix/combined_full.beagle.gz \
        -K $K -P 4 -seed 1 -minMaf 0.01 \
        -outfiles "$RESULTS_V1/K${K}" \
        > "$RESULTS_V1/K${K}.runlog" 2>&1
done

# v2: K=2..11
for K in $(seq 2 11); do
    log "NGSadmix v2: K=$K"
    "$NGSADMIX" -likes out_global/angsd/ngsadmix/combined_v2_full.beagle.gz \
        -K $K -P 4 -seed 1 -minMaf 0.01 \
        -outfiles "$RESULTS_V2/K${K}" \
        > "$RESULTS_V2/K${K}.runlog" 2>&1
done

# ── Step 4: Run PCAngsd ──────────────────────────────────────────────────────
log "=== Running PCAngsd ==="
log "NOTE: each pool = 1 diploid pseudo-individual in the Beagle (see warning above)"

PCANGSD_DIR="out_global/angsd/pcangsd"

pcangsd -b out_global/angsd/ngsadmix/combined_full.beagle.gz \
    -e 5 -t 4 --maf 0.01 \
    -o "$PCANGSD_DIR/pca_combined_full" \
    > "$PCANGSD_DIR/pca_combined_full.log" 2>&1
log "PCAngsd v1 complete"

pcangsd -b out_global/angsd/ngsadmix/combined_v2_full.beagle.gz \
    -e 5 -t 4 --maf 0.01 \
    -o "$PCANGSD_DIR/pca_combined_v2_full" \
    > "$PCANGSD_DIR/pca_combined_v2_full.log" 2>&1
log "PCAngsd v2 complete"

# ── Step 5: Run Evanno (multiple seeds) ──────────────────────────────────────
log "=== Running Evanno replicates ==="

EVANNO_V1="out_global/angsd/ngsadmix/results_evanno_v1_full"
EVANNO_V2="out_global/angsd/ngsadmix/results_evanno_v2_full"
mkdir -p "$EVANNO_V1" "$EVANNO_V2"

# v1: K=2..5, 4 seeds
for K in 2 3 4 5; do
    for SEED in 2 3 4 5; do
        log "Evanno v1: K=$K seed=$SEED"
        "$NGSADMIX" -likes out_global/angsd/ngsadmix/combined_full.beagle.gz \
            -K $K -P 3 -seed $SEED -minMaf 0.01 \
            -outfiles "$EVANNO_V1/K${K}_seed${SEED}" \
            > "$EVANNO_V1/K${K}_seed${SEED}.runlog" 2>&1
    done
done

# v2: K=2..11, 4 seeds
for K in $(seq 2 11); do
    for SEED in 2 3 4 5; do
        log "Evanno v2: K=$K seed=$SEED"
        "$NGSADMIX" -likes out_global/angsd/ngsadmix/combined_v2_full.beagle.gz \
            -K $K -P 3 -seed $SEED -minMaf 0.01 \
            -outfiles "$EVANNO_V2/K${K}_seed${SEED}" \
            > "$EVANNO_V2/K${K}_seed${SEED}.runlog" 2>&1
    done
done

# ── Step 6: Evanno deltaK calculation ────────────────────────────────────────
log "=== Calculating deltaK ==="
Rscript bin/10_evanno_deltaK.R "$EVANNO_V1" out_global/angsd/ngsadmix/evanno_v1_full
Rscript bin/10_evanno_deltaK.R "$EVANNO_V2" out_global/angsd/ngsadmix/evanno_v2_full

log "=== ALL DONE — Full genome panels complete ==="
