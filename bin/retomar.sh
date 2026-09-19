#!/usr/bin/env bash
set -euo pipefail

PROJ="$HOME/Documentos/Ancestria"
REF="$PROJ/data/ref/hg19.fa"
OUT="$PROJ/work/gatk_genome"
CHAIN="$PROJ/data/liftover/hg19ToHg38.over.chain.gz"
HG38="$PROJ/data/ref/hg38.fa"
mkdir -p "$PROJ/work/tmp"

run_sample() {
  local SAMPLE="$1" BAM="$2" PLOIDY="$3"
  echo "=== $SAMPLE (ploidía $PLOIDY) ==="

  pids=()
  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    vcf="$OUT/${SAMPLE}_chr${chr}_raw.vcf.gz"
    if [[ -s "$vcf" ]]; then
      echo "  [=] $SAMPLE chr${chr} ya existe"
      continue
    fi
    echo "  [..] $SAMPLE chr${chr}..."
    (
      gatk HaplotypeCaller \
        -R "$REF" -I "$BAM" -O "$vcf" -L "$chr" \
        --sample-ploidy "$PLOIDY" \
        --native-pair-hmm-threads 1 \
        --tmp-dir "$PROJ/work/tmp" -ERC NONE \
        > "$OUT/${SAMPLE}_chr${chr}.log" 2>&1 \
      && echo "  [✓] $SAMPLE chr${chr} OK" \
      || echo "  [✗] $SAMPLE chr${chr} FALLÓ"
    ) &
    pids+=($!)
    if [[ ${#pids[@]} -ge 11 ]]; then
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
    fi
  done
  for pid in "${pids[@]}"; do wait "$pid"; done

  echo "  [...] Merge $SAMPLE..."
  vcf_list=()
  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    vcf_list+=("$OUT/${SAMPLE}_chr${chr}_raw.vcf.gz")
  done
  bcftools concat -a -D -O z -o "$OUT/${SAMPLE}_merged_raw.vcf.gz" "${vcf_list[@]}"
  bcftools index -t "$OUT/${SAMPLE}_merged_raw.vcf.gz"

  echo "  [...] Filtrado $SAMPLE..."
  bcftools view -m2 -M2 -v snps \
    -i 'QUAL>=30 && INFO/DP>=20 && INFO/AF>=0.01 && REF!="N" && ALT!="."' \
    -O z -o "$OUT/${SAMPLE}_filtered.vcf.gz" \
    "$OUT/${SAMPLE}_merged_raw.vcf.gz"
  bcftools index -t "$OUT/${SAMPLE}_filtered.vcf.gz"

  echo "  [...] Liftover $SAMPLE..."
  bcftools view "$OUT/${SAMPLE}_filtered.vcf.gz" | \
    iconv -f utf-8 -t utf-8 -c > "$OUT/${SAMPLE}_clean.vcf"
  CrossMap vcf "$CHAIN" "$OUT/${SAMPLE}_clean.vcf" "$HG38" \
    "$OUT/${SAMPLE}_hg38_unsorted.vcf"
  bcftools sort -O z -o "$OUT/${SAMPLE}_hg38.vcf.gz" \
    "$OUT/${SAMPLE}_hg38_unsorted.vcf"
  bcftools index -t "$OUT/${SAMPLE}_hg38.vcf.gz"
  rm -f "$OUT/${SAMPLE}_clean.vcf" "$OUT/${SAMPLE}_hg38_unsorted.vcf"

  N=$(bcftools stats "$OUT/${SAMPLE}_hg38.vcf.gz" | grep "number of SNPs:" | cut -f4)
  echo "  [✓] $SAMPLE completo: $N SNPs en hg38"
}

# Lanzar ambos en paralelo
run_sample POOL2 "$PROJ/data/bam/51-100_samples_UDB-103_482264.merged.bam" 100 \
  > "$PROJ/logs/POOL2_gatk.log" 2>&1 &

run_sample HOSPITAL "$PROJ/data/bam/HOSPITAL.bam" 400 \
  > "$PROJ/logs/HOSPITAL_gatk.log" 2>&1 &

wait
echo "=== TODO COMPLETADO ==="
