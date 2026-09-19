#!/usr/bin/env bash
# 18_build_X_panel_and_run.sh — builds a de novo NGSadmix/PCAngsd panel for
# chromosome X and runs NGSadmix + PCAngsd on it.
#
# CONTEXT: the HapMap3 panel used for autosomes (bin/06_build_5pop_panel.sh,
# bin/11_full_genome_panels.sh) has 0 positions on X. This script defines an
# alternative panel from scratch:
#   1. Sites: intersection of ANGSD's independent major/minor calls across
#      the 3 pools (requires bin/05_angsd_analysis.sh --chr X to have run
#      first).
#   2. Reference: remote streaming of 1000G Phase 3 chrX (no full VCF
#      download -- a single sequential pass + local filtering, much faster
#      than querying region by region with `bcftools -R`).
#   3. Pools + reference Beagle, merge, NGSadmix K=5, PCAngsd.
#
# WARNING (see docs/NGSADMIX_X_STATUS.md): this script's
# result does NOT pass validation (the 3 pools come out nearly identical,
# contradicting chromosome X's own FST computed from the same data, and
# PCAngsd does not converge). The de novo panel built here has only ~4,350
# SNPs, far below the ~245,000 of the autosomal panel, and apparently is not
# enough to resolve population structure at this pool-seq's depth. The
# script is kept for reproducibility and as a starting point for a retry
# with better marker selection (see status.md).
#
# Requires: bin/05_angsd_analysis.sh --chr X already run (produces the
# *_maf.mafs.gz used to define sites), individuals_full.tsv from the
# autosomal run (bin/12_rebuild_and_run.sh), bcftools, angsd, NGSadmix,
# pcangsd.
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ_DIR"

WORK="work/angsd_X_panel"
OUT="out_global/angsd_X/ngsadmix"
ANGSD="$PROJ_DIR/tools/angsd_src/angsd"
NGSADMIX="$PROJ_DIR/tools/angsd_src/misc/NGSadmix"
PCANGSD="${PCANGSD_BIN:-pcangsd}"  # set PCANGSD_BIN if not on PATH
REF_FASTA="data/ref/human_g1k_v37_decoy.fasta"
XVCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"

mkdir -p "$WORK/ref_extract_final" "$OUT"
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Step 1: de novo sites from the 3 pools' MAF on X ────────────────────────
log "Step 1/6: defining sites from the pools' MAF"
for p in POOL1 POOL2 HOSPITAL; do
    f="out_global/angsd_X/per_pool/${p}_maf.mafs.gz"
    [[ -f "$f" ]] || { echo "ERROR: missing $f. Run first: bash bin/05_angsd_analysis.sh --ref $REF_FASTA --chr X" >&2; exit 1; }
done
Rscript -e '
suppressPackageStartupMessages(library(data.table))
p1 <- fread("out_global/angsd_X/per_pool/POOL1_maf.mafs.gz")
p2 <- fread("out_global/angsd_X/per_pool/POOL2_maf.mafs.gz")
p3 <- fread("out_global/angsd_X/per_pool/HOSPITAL_maf.mafs.gz")
setnames(p1, c("chromo","position","major","minor"), c("chr","pos","maj1","min1"))
setnames(p2, c("chromo","position","major","minor"), c("chr","pos","maj2","min2"))
setnames(p3, c("chromo","position","major","minor"), c("chr","pos","maj3","min3"))
m <- merge(p1[,.(chr,pos,maj1,min1)], p2[,.(chr,pos,maj2,min2)], by=c("chr","pos"))
m <- merge(m, p3[,.(chr,pos,maj3,min3)], by=c("chr","pos"))
consistent <- m[maj1==maj2 & maj2==maj3 & min1==min2 & min2==min3]
out <- consistent[, .(chr, pos, major=maj1, minor=min1)]
setorder(out, chr, pos)
fwrite(out, "work/angsd_X_panel/sites_X.txt", sep="\t", col.names=FALSE)
cat(sprintf("Sites consistent across the 3 pools: %d\n", nrow(out)))
'
"$ANGSD" sites index "$WORK/sites_X.txt"

# ── Step 2: pools' Beagle at those sites ─────────────────────────────────────
log "Step 2/6: pools' genotype likelihoods (doGlf 2)"
printf '%s\n' \
    bam/01-50_samples_UDB-102_482263.merged.bam \
    bam/51-100_samples_UDB-103_482264.merged.bam \
    bam/HOSPITAL.bam > "$WORK/pools.bamlist"
"$ANGSD" -bam "$WORK/pools.bamlist" -r X: \
    -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 20 \
    -ref "$REF_FASTA" -sites "$WORK/sites_X.txt" -checkBamHeaders 0 \
    -GL 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -nThreads 4 \
    -out "$OUT/pools_X"

# ── Step 3: reference genotypes via remote streaming (no full download) ────
log "Step 3/6: reference genotypes (remote streaming of 1000G chrX)"
tail -n +2 out_global/angsd/ngsadmix/individuals_full.tsv | grep -v "^POOL\|^HOSPITAL" | cut -f1 > "$WORK/ref_samples.txt"
awk '{print $1"\t"$2}' "$WORK/sites_X.txt" | sort -k2,2n > "$WORK/wanted_pos.txt"
# NOTE: a single sequential pass over all of chrX + local filtering is MUCH
# faster than `bcftools -R` with thousands of small regions (each one
# triggers a separate remote query). Do not use -R with many regions over a
# remote URL.
bcftools view -r X -S "$WORK/ref_samples.txt" --force-samples -v snps -m2 -M2 "$XVCF_URL" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  | awk -v wanted="$WORK/wanted_pos.txt" \
        'BEGIN{while((getline line < wanted)>0){split(line,a,"\t"); want[a[2]]=1}} ($2 in want)' \
  > "$WORK/ref_extract_final/chrX_gt.txt"

# ── Step 4: reference Beagle + merge ─────────────────────────────────────────
log "Step 4/6: reference Beagle and merge with pools"
Rscript "$PROJ_DIR/bin/lib_build_ref_beagle_panel.R" \
    "$WORK/sites_X.txt" "$WORK/ref_extract_final" "$WORK/ref_samples.txt" \
    "$OUT/refX.beagle" "X (100 individuals)" X

Rscript -e '
suppressPackageStartupMessages(library(data.table))
POOL_NAMES <- c("POOL1","POOL2","HOSPITAL")
pools <- fread("out_global/angsd_X/ngsadmix/pools_X.beagle.gz", header=TRUE)
ind_idx <- rep(seq_along(POOL_NAMES), each=3)
setnames(pools, 4:ncol(pools), paste0(POOL_NAMES[ind_idx], c("_p1","_p2","_p3")))
ref <- fread("out_global/angsd_X/ngsadmix/refX.beagle", header=FALSE)
setnames(ref, 1:3, c("marker","allele1","allele2"))
m <- merge(pools, ref, by=c("marker","allele1","allele2"))
fwrite(m, "out_global/angsd_X/ngsadmix/combined_X.beagle.gz", sep="\t")
samples <- readLines("work/angsd_X_panel/ref_samples.txt")
panel <- fread("data/panel/integrated_call_samples_v3.20130502.ALL.panel", fill=TRUE, skip=1)
setnames(panel, 1:4, c("sample","pop","super_pop","gender"))
pop_map <- panel[match(samples, sample), .(sample, pop, super_pop)]
ind_table <- rbind(data.table(sample=POOL_NAMES, pop=POOL_NAMES, super_pop="POOL"), pop_map)
fwrite(ind_table, "out_global/angsd_X/ngsadmix/individuals_X.tsv", sep="\t")
cat(sprintf("Combined: %d sites, %d individuals\n", nrow(m), (ncol(m)-3)/3))
'

# ── Step 5: NGSadmix K=5 ──────────────────────────────────────────────────────
log "Step 5/6: NGSadmix K=5"
mkdir -p "$OUT/results"
"$NGSADMIX" -likes "$OUT/combined_X.beagle.gz" -K 5 -P 4 -seed 1 -minMaf 0.01 \
    -outfiles "$OUT/results/K5_X" > "$OUT/results/K5_X.runlog" 2>&1

# ── Step 6: PCAngsd ───────────────────────────────────────────────────────────
log "Step 6/6: PCAngsd"
mkdir -p out_global/angsd_X/pcangsd
"$PCANGSD" -b "$OUT/combined_X.beagle.gz" -e 5 -t 4 --maf 0.01 \
    -o out_global/angsd_X/pcangsd/pca_X > out_global/angsd_X/pcangsd/pca_X.log 2>&1

log "Done. SEE docs/NGSADMIX_X_STATUS.md before using these results."
