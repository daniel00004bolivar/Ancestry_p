#!/usr/bin/env bash
# 19a_prepare_y_markers_and_counts.sh — prepares the inputs for
# bin/19_y_haplogroup_mixture.R: downloads the 1000G Phase 3 Y VCF, selects
# SNPs with strong cross-superpopulation differentiation, and runs ANGSD
# ACGT counts on each pool at those sites.
#
# WARNING (see docs/Y_HAPLOGROUP_STATUS.md): this pipeline's
# final result (bin/19_y_haplogroup_mixture.R) does NOT pass validation --
# too few sites with sufficient depth, strongly correlated with each other
# (Y does not recombine), and the result is unstable and implausible across
# pools. Kept for reproducibility and as a starting point if more Y-specific
# sequencing depth becomes available.
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ_DIR"

WORK="work/y_panel"
ANGSD="$PROJ_DIR/tools/angsd_src/angsd"
REF_FASTA="data/ref/human_g1k_v37_decoy.fasta"
YVCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
SPREAD_MIN="${SPREAD_MIN:-0.6}"  # cross-superpopulation differentiation threshold

mkdir -p "$WORK" out_global/angsd_Y/per_pool
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Step 1: download the Y VCF (small, ~5.4MB, no streaming needed) ────────
log "Step 1/4: downloading the 1000G chrY VCF"
curl -sL -o "$WORK/chrY.vcf.gz" "$YVCF_URL"
curl -sL -o "$WORK/chrY.vcf.gz.tbi" "${YVCF_URL}.tbi"

# ── Step 2: extract per-superpopulation frequencies from the VCF's INFO ────
log "Step 2/4: extracting per-superpopulation frequencies"
bcftools view -H -v snps -m2 -M2 "$WORK/chrY.vcf.gz" | awk 'BEGIN{OFS="\t"}
{
    chr=$1; pos=$2; ref=$4; alt=$5; info=$8
    if (length(ref)!=1 || length(alt)!=1) next
    amr="."; afr="."; eas="."; eur="."; sas="."
    n=split(info, f, ";")
    for(i=1;i<=n;i++) {
        if (f[i] ~ /^AMR_AF=/) amr=substr(f[i],8)
        if (f[i] ~ /^AFR_AF=/) afr=substr(f[i],8)
        if (f[i] ~ /^EAS_AF=/) eas=substr(f[i],8)
        if (f[i] ~ /^EUR_AF=/) eur=substr(f[i],8)
        if (f[i] ~ /^SAS_AF=/) sas=substr(f[i],8)
    }
    if (amr=="."||afr=="."||eas=="."||eur=="."||sas==".") next
    print chr, pos, ref, alt, afr, amr, eas, eur, sas
}' > "$WORK/y_pop_freqs.txt"

# ── Step 3: select the most differentiated markers and build the sites file ─
log "Step 3/4: selecting markers with spread > $SPREAD_MIN"
awk -v thr="$SPREAD_MIN" 'BEGIN{OFS="\t"}
{
    afr=$5; amr=$6; eas=$7; eur=$8; sas=$9
    max=afr; if(amr>max)max=amr; if(eas>max)max=eas; if(eur>max)max=eur; if(sas>max)max=sas
    min=afr; if(amr<min)min=amr; if(eas<min)min=eas; if(eur<min)min=eur; if(sas<min)min=sas
    spread=max-min
    print $0, spread
}' "$WORK/y_pop_freqs.txt" | sort -k10,10gr > "$WORK/y_freqs_sorted.txt"
awk -v thr="$SPREAD_MIN" '$10>thr {print $1"\t"$2"\t"$3"\t"$4}' "$WORK/y_freqs_sorted.txt" > "$WORK/y_markers.txt"
wc -l "$WORK/y_markers.txt"
"$ANGSD" sites index "$WORK/y_markers.txt"

# ── Step 4: ANGSD ACGT counts per pool ───────────────────────────────────────
log "Step 4/4: per-pool ACGT counts"
declare -A BAMS=(
    [POOL1]="bam/01-50_samples_UDB-102_482263.merged.bam"
    [POOL2]="bam/51-100_samples_UDB-103_482264.merged.bam"
    [HOSPITAL]="bam/HOSPITAL.bam"
)
for name in POOL1 POOL2 HOSPITAL; do
    "$ANGSD" -i "${BAMS[$name]}" -r Y: \
        -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 20 \
        -ref "$REF_FASTA" -sites "$WORK/y_markers.txt" \
        -doCounts 1 -dumpCounts 3 -nThreads 4 \
        -out "out_global/angsd_Y/per_pool/${name}_counts"
done

log "Done. Next: Rscript bin/19_y_haplogroup_mixture.R"
log "SEE docs/Y_HAPLOGROUP_STATUS.md before using these results."
