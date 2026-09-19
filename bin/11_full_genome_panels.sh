#!/usr/bin/env bash
set -euo pipefail

PROJ_DIR="/Users/danielbolivar/Downloads/Ancestria"
cd "$PROJ_DIR"

SAMPLES_V1="/tmp/ref_sample_order.txt"
SAMPLES_V2="/tmp/ref_sample_order_v2.txt"
SITES="/tmp/angsd_sites_hapmap3.txt"
PANEL_FILE="data/panel/integrated_call_samples_v3.20130502.ALL.panel"

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_DIR="/tmp/1kgp_vcfs"
EXTRACT_V1="/tmp/ref_extract"
EXTRACT_V2="/tmp/ref_extract_v2"

mkdir -p "$VCF_DIR" "$EXTRACT_V1" "$EXTRACT_V2"

ALL_CHRS=$(seq 1 22)
EXISTING_CHRS="1 6 10 13 17 20 22"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Step 1: Download missing VCFs ────────────────────────────────────────────
for CHR in $ALL_CHRS; do
    VCF_FILE="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    VCF_PATH="$VCF_DIR/$VCF_FILE"
    TBI_PATH="${VCF_PATH}.tbi"

    if [[ -f "$VCF_PATH" ]]; then
        log "chr$CHR: VCF already downloaded, skipping"
        continue
    fi

    log "chr$CHR: Downloading VCF..."
    curl -L -o "$VCF_PATH" "${BASE_URL}/${VCF_FILE}" 2>/dev/null &
    CURL_PID=$!

    if [[ ! -f "$TBI_PATH" ]]; then
        # Check if .tbi is in PROJ_DIR
        LOCAL_TBI="$PROJ_DIR/$VCF_FILE.tbi"
        if [[ -f "$LOCAL_TBI" ]]; then
            cp "$LOCAL_TBI" "$TBI_PATH"
            log "chr$CHR: TBI copied from project dir"
        else
            log "chr$CHR: Downloading TBI..."
            curl -L -o "$TBI_PATH" "${BASE_URL}/${VCF_FILE}.tbi" 2>/dev/null
        fi
    fi

    wait $CURL_PID
    log "chr$CHR: Download complete ($(du -h "$VCF_PATH" | cut -f1))"
done

# ── Step 2: Build regions files for all chromosomes ──────────────────────────
log "Building regions files for all chromosomes..."
for CHR in $ALL_CHRS; do
    REGIONS_V1="$EXTRACT_V1/regions_chr${CHR}.txt"
    REGIONS_V2="$EXTRACT_V2/regions_chr${CHR}.txt"
    if [[ ! -f "$REGIONS_V1" ]]; then
        awk -v c="$CHR" '$1 == c {print $1 "\t" $2}' "$SITES" > "$REGIONS_V1"
    fi
    if [[ ! -f "$REGIONS_V2" ]]; then
        cp "$REGIONS_V1" "$REGIONS_V2"
    fi
done

# ── Step 3: Extract genotypes for all chromosomes ────────────────────────────
for CHR in $ALL_CHRS; do
    VCF_FILE="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    VCF_PATH="$VCF_DIR/$VCF_FILE"
    TBI_PATH="${VCF_PATH}.tbi"

    # Ensure TBI exists
    if [[ ! -f "$TBI_PATH" ]]; then
        LOCAL_TBI="$PROJ_DIR/$VCF_FILE.tbi"
        if [[ -f "$LOCAL_TBI" ]]; then
            cp "$LOCAL_TBI" "$TBI_PATH"
        else
            log "chr$CHR: Downloading TBI..."
            curl -L -o "$TBI_PATH" "${BASE_URL}/${VCF_FILE}.tbi" 2>/dev/null
        fi
    fi

    # Panel v1 (100 individuals)
    GT_V1="$EXTRACT_V1/chr${CHR}_gt.txt"
    if [[ ! -f "$GT_V1" ]] || [[ ! -s "$GT_V1" ]]; then
        REGIONS_V1="$EXTRACT_V1/regions_chr${CHR}.txt"
        N_SITES=$(wc -l < "$REGIONS_V1")
        if [[ "$N_SITES" -gt 0 ]]; then
            log "chr$CHR: Extracting genotypes for 100 individuals ($N_SITES sites)..."
            bcftools view -R "$REGIONS_V1" -S "$SAMPLES_V1" --force-samples \
                -v snps -m2 -M2 "$VCF_PATH" \
              | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
              > "$GT_V1" 2>"$EXTRACT_V1/chr${CHR}_err.log"
            log "chr$CHR: v1 extracted $(wc -l < "$GT_V1") sites"
        fi
    else
        log "chr$CHR: v1 genotypes already exist ($(wc -l < "$GT_V1") sites)"
    fi

    # Panel v2 (110 individuals)
    GT_V2="$EXTRACT_V2/chr${CHR}_gt.txt"
    if [[ ! -f "$GT_V2" ]] || [[ ! -s "$GT_V2" ]]; then
        REGIONS_V2="$EXTRACT_V2/regions_chr${CHR}.txt"
        N_SITES=$(wc -l < "$REGIONS_V2")
        if [[ "$N_SITES" -gt 0 ]]; then
            log "chr$CHR: Extracting genotypes for 110 individuals ($N_SITES sites)..."
            bcftools view -R "$REGIONS_V2" -S "$SAMPLES_V2" --force-samples \
                -v snps -m2 -M2 "$VCF_PATH" \
              | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
              > "$GT_V2" 2>"$EXTRACT_V2/chr${CHR}_err.log"
            log "chr$CHR: v2 extracted $(wc -l < "$GT_V2") sites"
        fi
    else
        log "chr$CHR: v2 genotypes already exist ($(wc -l < "$GT_V2") sites)"
    fi
done

log "=== All genotype extraction complete ==="
for CHR in $ALL_CHRS; do
    V1=$(wc -l < "$EXTRACT_V1/chr${CHR}_gt.txt" 2>/dev/null || echo 0)
    V2=$(wc -l < "$EXTRACT_V2/chr${CHR}_gt.txt" 2>/dev/null || echo 0)
    echo "  chr$CHR: v1=$V1 sites, v2=$V2 sites"
done
