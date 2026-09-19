#!/usr/bin/env bash
# 05_angsd_analysis.sh — ANGSD analysis for pool-seq ancestry
#
# Analyzes 3 pool-seq BAMs (POOL1, POOL2, HOSPITAL) with ANGSD.
# No reference genome required for the basic analyses (1-3).
# Analysis 4 (SFS/FST via realSFS) requires --ref.
#
# Analyses:
#   1. Per-pool allele frequencies (doMaf)
#   2. Per-pool allele counts (doCounts)
#   3. IBS matrix across the 3 pools
#   4. [With --ref] SAF + 1D SFS + pairwise FST via realSFS
#
# Pool-seq note: ANGSD assumes each BAM is one diploid individual. For
# pool-seq, GL is biased, but counts (doCounts) and frequencies (doMaf) are
# a reasonable approximation given sufficient depth. Counts-based FST
# (Hudson) is the most robust estimator for pools.
#
# EXTENDED METHODOLOGICAL WARNING: this is not just a GL bias in this
# script -- it propagates downstream to NGSadmix/PCAngsd (bin/08*, bin/09*),
# which represent each pool as ONE diploid pseudo-individual with no way to
# indicate its real size (~50 individuals). See bin/README_VIGENTE.md and
# bin/16_pool_replicates.sh (uncertainty quantification via block jackknife).
#
# Usage:
#   bash bin/05_angsd_analysis.sh [options]
#
# Options:
#   --ref    /path/to/human_g1k_v37_decoy.fasta   (enables SFS and FST via realSFS)
#   --chr    1-22  (chromosomes; default: 1-22)
#   --region 21,22  (fast test region; e.g. "21,22")
#   --threads N     (default: 4)
#   --minQ   N      (minimum base quality; default: 20)
#   --minMQ  N      (minimum mapping quality; default: 30)
#   --out    dir    (output directory; default: out_global/angsd)
#   --dry-run       (print commands without running)
set -euo pipefail

# ── Rutas ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ANGSD="$PROJ_DIR/tools/angsd_src/angsd"
REALSFS="$PROJ_DIR/tools/angsd_src/misc/realSFS"

BAM_POOL1="$PROJ_DIR/bam/01-50_samples_UDB-102_482263.merged.bam"
BAM_POOL2="$PROJ_DIR/bam/51-100_samples_UDB-103_482264.merged.bam"
BAM_HOSP="$PROJ_DIR/bam/HOSPITAL.bam"

# ── Valores por defecto ────────────────────────────────────────────────────────
REF=""
THREADS=4
MIN_Q=20
MIN_MQ=30
OUT_DIR="$PROJ_DIR/out_global/angsd"
CHR_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
DRY_RUN=false

# ── Parse args ─────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --ref)     REF="$2";     shift 2 ;;
        --chr)     CHR_LIST="$2"; shift 2 ;;
        --region)  CHR_LIST="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --minQ)    MIN_Q="$2";   shift 2 ;;
        --minMQ)   MIN_MQ="$2";  shift 2 ;;
        --out)     OUT_DIR="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── Validation ───────────────────────────────────────────────────────────────
for bam in "$BAM_POOL1" "$BAM_POOL2" "$BAM_HOSP"; do
    [[ -f "$bam" ]] || { echo "ERROR: BAM not found: $bam" >&2; exit 1; }
    [[ -f "${bam}.bai" ]] || { echo "ERROR: BAM index not found: ${bam}.bai" >&2; exit 1; }
done
[[ -x "$ANGSD" ]]   || { echo "ERROR: angsd not found at $ANGSD" >&2; exit 1; }
[[ -x "$REALSFS" ]] || { echo "ERROR: realSFS not found at $REALSFS" >&2; exit 1; }
if [[ -n "$REF" ]] && [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference not found: $REF" >&2; exit 1
fi

mkdir -p "$OUT_DIR"/{per_pool,ibs,fst,sfs}

# ── Helpers ────────────────────────────────────────────────────────────────────
run() {
    if $DRY_RUN; then
        echo "[DRY-RUN] $*"
    else
        echo "[CMD] $*"
        "$@"
    fi
}

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# Converts "1,2,22" or "1-22" or "21,22" into a region list for -r
build_rf_list() {
    local spec="$1"
    local tmpfile="$OUT_DIR/regions.txt"
    # Range (has a dash): expand it
    if [[ "$spec" == *-* && "$spec" != *,* ]]; then
        local s="${spec%-*}" e="${spec#*-}"
        for ((i=s; i<=e; i++)); do echo "$i"; done > "$tmpfile"
    else
        # Comma-separated list
        echo "$spec" | tr ',' '\n' > "$tmpfile"
    fi
    echo "$tmpfile"
}

# ── Base filters (common to all analyses) ───────────────────────────────────
# -C 50 and -baq 1 require a reference; omitted if --ref was not given.
# FIXED BUG (audit): when --ref is given, BASE_FILTERS adds -C/-baq, but none
# of blocks 1-3 (MAF, counts, IBS) passed -ref to the angsd command -- only
# block 4 (SAF) did. angsd fails with "Must also supply -ref for adjusting
# the mapping quality" and produces empty output WITHOUT stopping the script
# (the error does not propagate as a failure under set -e in this context).
# REF_FLAG is now also added to blocks 1-3 when applicable.
if [[ -n "$REF" ]]; then
    BASE_FILTERS="-remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ $MIN_MQ -minQ $MIN_Q"
    REF_FLAG="-ref $REF"
else
    BASE_FILTERS="-remove_bads 1 -only_proper_pairs 1 -minMapQ $MIN_MQ -minQ $MIN_Q"
    REF_FLAG=""
fi

# ── BAM lists ─────────────────────────────────────────────────────────────────
BAM_LIST_ALL="$OUT_DIR/all_pools.bamlist"
printf '%s\n' "$BAM_POOL1" "$BAM_POOL2" "$BAM_HOSP" > "$BAM_LIST_ALL"

RF_FILE=$(build_rf_list "$CHR_LIST")

log "=== ANGSD Pool-Seq Ancestry Analysis ==="
log "Pools: POOL1 | POOL2 | HOSPITAL"
log "Chromosomes: $CHR_LIST"
log "Reference: ${REF:-'not provided (no-ref mode)'}"
log "Output: $OUT_DIR"
log "Threads: $THREADS"

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Per-pool allele frequencies (no reference)
# doMajorMinor 1 = infer major/minor from the data
# doMaf 1        = compute alternate allele frequency
# GL 1           = SAMtools-style GL model
# SNP_pval 1e-6  = filter monomorphic sites
# ══════════════════════════════════════════════════════════════════════════════
log "--- [1/4] Estimating per-pool allele frequencies ---"
# Pool-seq note: -SNP_pval needs N>1 ind for a chi2 test; we use -minMaf instead.
# -doMajorMinor 1 infers the major/minor allele from GL (SAMtools model, -GL 1).

bam_for_sample() {
    case "$1" in
        POOL1)    echo "$BAM_POOL1" ;;
        POOL2)    echo "$BAM_POOL2" ;;
        HOSPITAL) echo "$BAM_HOSP"  ;;
    esac
}

for SAMPLE in POOL1 POOL2 HOSPITAL; do
    BAM="$(bam_for_sample "$SAMPLE")"
    OUT="$OUT_DIR/per_pool/${SAMPLE}_maf"
    log "  Processing $SAMPLE → $OUT"
    run "$ANGSD" \
        -i "$BAM" \
        -rf "$RF_FILE" \
        $BASE_FILTERS \
        $REF_FLAG \
        -GL 1 \
        -doMajorMinor 1 \
        -doMaf 1 \
        -minMaf 0.01 \
        -nThreads "$THREADS" \
        -out "$OUT"
done

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Per-pool allele counts (A, C, G, T per site)
# doCounts 1     = count reads per allele
# dumpCounts 2   = dump ref/alt counts (2 cols: ref_n, alt_n)
# doMajorMinor 4 = infer major/minor from counts (not GL), appropriate for pools
# These counts underlie the Hudson FST estimator (more robust than GL for pools)
# ══════════════════════════════════════════════════════════════════════════════
log "--- [2/4] Dumping per-pool allele counts ---"

for SAMPLE in POOL1 POOL2 HOSPITAL; do
    BAM="$(bam_for_sample "$SAMPLE")"
    OUT="$OUT_DIR/per_pool/${SAMPLE}_counts"
    log "  Processing $SAMPLE → $OUT"
    run "$ANGSD" \
        -i "$BAM" \
        -rf "$RF_FILE" \
        $BASE_FILTERS \
        $REF_FLAG \
        -doCounts 1 \
        -dumpCounts 3 \
        -setMinDepth 5 \
        -nThreads "$THREADS" \
        -out "$OUT"
done
# Per-pool outputs: ${SAMPLE}_counts.pos.gz (chr, pos, totalDepth)
#                   ${SAMPLE}_counts.counts.gz (totA, totC, totG, totT)

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: IBS matrix across the 3 pools (no reference)
# doIBS 1           = pairwise identity-by-state (3×3 matrix)
# doCounts 1        = required for doIBS
# checkBamHeaders 0 = ignore contig differences between POOL1/2 and HOSPITAL
#                     (HOSPITAL has extra NC_007605 and hs37d5 contigs)
# NOTE: for pool-seq, each pool is treated as a "pseudo-individual".
#        IBS distance between pools reflects population differentiation.
# ══════════════════════════════════════════════════════════════════════════════
log "--- [3/4] Computing IBS matrix across pools ---"
OUT_IBS="$OUT_DIR/ibs/all_pools"
run "$ANGSD" \
    -bam "$BAM_LIST_ALL" \
    -rf "$RF_FILE" \
    $BASE_FILTERS \
    $REF_FLAG \
    -doCounts 1 \
    -doIBS 1 \
    -doMajorMinor 4 \
    -doMaf 1 \
    -minMaf 0.01 \
    -checkBamHeaders 0 \
    -nThreads "$THREADS" \
    -out "$OUT_IBS"

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: 1D SFS + pairwise FST via realSFS (REQUIRES REFERENCE)
# doSaf 1   = site allele frequency likelihoods (needs anc/ref)
# ══════════════════════════════════════════════════════════════════════════════
if [[ -n "$REF" ]]; then
    log "--- [4/4] 1D SFS and pairwise FST via realSFS (with reference) ---"

    for SAMPLE in POOL1 POOL2 HOSPITAL; do
        BAM="$(bam_for_sample "$SAMPLE")"
        OUT="$OUT_DIR/sfs/${SAMPLE}_saf"
        log "  SAF for $SAMPLE → $OUT"
        run "$ANGSD" \
            -i "$BAM" \
            -rf "$RF_FILE" \
            $BASE_FILTERS \
            -GL 1 \
            -doSaf 1 \
            -anc "$REF" \
            -ref "$REF" \
            -nThreads "$THREADS" \
            -out "$OUT"

        log "  1D SFS for $SAMPLE"
        run "$REALSFS" "${OUT}.saf.idx" -P "$THREADS" > "${OUT}.sfs"
    done

    # Pairwise FST
    log "  Pairwise FST via realSFS"
    PAIRS=("POOL1 POOL2" "POOL1 HOSPITAL" "POOL2 HOSPITAL")
    for pair in "${PAIRS[@]}"; do
        A="${pair% *}"; B="${pair#* }"
        SAF_A="$OUT_DIR/sfs/${A}_saf.saf.idx"
        SAF_B="$OUT_DIR/sfs/${B}_saf.saf.idx"
        ML_OUT="$OUT_DIR/fst/${A}_${B}.ml"
        FST_OUT="$OUT_DIR/fst/${A}_${B}"
        log "    FST: $A vs $B"
        run "$REALSFS" "$SAF_A" "$SAF_B" -P "$THREADS" > "$ML_OUT"
        run "$REALSFS" fst index "$SAF_A" "$SAF_B" -sfs "$ML_OUT" -fstout "$FST_OUT" -P "$THREADS"
        run "$REALSFS" fst stats "${FST_OUT}.fst.idx"
    done
else
    log "--- [4/4] SKIPPED: SFS/FST via realSFS (pass --ref to enable it) ---"
    log "    To download the b37 reference: see the README at the end of this script"
fi

# ══════════════════════════════════════════════════════════════════════════════
# POST-PROCESSING: Hudson FST from counts (no reference needed)
# Produces an FST table for the 3 pool pairs using allele counts.
# Hudson et al. 1992 estimator: Fst = 1 - (Hs/Ht)
# Only runs if not --dry-run and if the counts files exist.
# ══════════════════════════════════════════════════════════════════════════════
if ! $DRY_RUN; then
    log "--- Post-proc: Hudson FST from counts (no reference) ---"
    Rscript "$SCRIPT_DIR/05a_fst_from_counts.R" \
        --counts_dir "$OUT_DIR/per_pool" \
        --out "$OUT_DIR/fst/hudson_fst_summary.tsv" \
        2>&1 | sed 's/^/    /' || log "  (R script pending or post-proc error)"
fi

log "=== ANGSD analysis complete ==="
log "Results in: $OUT_DIR"
log ""
log "Key files:"
log "  Frequencies: $OUT_DIR/per_pool/*_maf.mafs.gz"
log "  Counts:      $OUT_DIR/per_pool/*_counts.counts.gz"
log "  IBS:         $OUT_DIR/ibs/all_pools.ibs.gz"
[[ -n "$REF" ]] && log "  FST realSFS: $OUT_DIR/fst/*.fst.idx"
log "  FST Hudson:  $OUT_DIR/fst/hudson_fst_summary.tsv"
log ""
log "To run a fast test on chr21-22:"
log "  bash bin/05_angsd_analysis.sh --region 21,22 --threads 4"
log ""
log "To run with a reference (full SFS + robust FST):"
log "  bash bin/05_angsd_analysis.sh --ref /path/to/human_g1k_v37_decoy.fasta"
log ""
log "b37 reference (if you don't have it):"
log "  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37_decoy.fasta.gz"
log "  bgzip -d human_g1k_v37_decoy.fasta.gz"
log "  samtools faidx human_g1k_v37_decoy.fasta"
