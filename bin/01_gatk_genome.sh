#!/usr/bin/env bash
# =============================================================================
# 01_gatk_genome.sh — v2
# Sin export -f. Cada job corre en subshell ( ) con su propio log.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

THREADS="${1:-20}"
REF_FASTA="$PROJ_DIR/data/ref/hg19.fa"
BAM_DIR="$PROJ_DIR/data/bam"
OUT_DIR="$WORK/gatk_genome"
CHAIN="$PROJ_DIR/data/liftover/hg19ToHg38.over.chain.gz"
HG38="$PROJ_DIR/data/ref/hg38.fa"
mkdir -p "$OUT_DIR" "$WORK/tmp"

declare -A BAMS=(
  [HOSPITAL]="$BAM_DIR/HOSPITAL.bam"
  [POOL1]="$BAM_DIR/01-50_samples_UDB-102_482263.merged.bam"
  [POOL2]="$BAM_DIR/51-100_samples_UDB-103_482264.merged.bam"
)
declare -A PLOIDY=([HOSPITAL]=400 [POOL1]=100 [POOL2]=100)

CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
PARALLEL_CHRS=10
CORES_PER_CHR=$(( THREADS / PARALLEL_CHRS ))
[[ "$CORES_PER_CHR" -lt 1 ]] && CORES_PER_CHR=1

echo "============================================================"
echo " GATK HaplotypeCaller — Genoma completo v2"
echo " Threads: $THREADS | Paralelo: $PARALLEL_CHRS chr | Cores/chr: $CORES_PER_CHR"
echo "============================================================"

for sample in POOL1 POOL2 HOSPITAL; do
  bam="${BAMS[$sample]}"
  ploidy="${PLOIDY[$sample]}"

  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo " $sample | ploidía $ploidy"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  pids=()
  for chr in "${CHRS[@]}"; do
    out_vcf="$OUT_DIR/${sample}_chr${chr}_raw.vcf.gz"

    if [[ -s "$out_vcf" ]]; then
      echo "  [=] $sample chr${chr} ya existe."
      continue
    fi

    echo "  [..] $sample chr${chr}..."

    (
      gatk HaplotypeCaller \
        -R "$REF_FASTA" \
        -I "$bam" \
        -O "$out_vcf" \
        -L "${chr}" \
        --sample-ploidy "$ploidy" \
        --native-pair-hmm-threads "$CORES_PER_CHR" \
        --tmp-dir "$WORK/tmp" \
        -ERC NONE \
        > "$OUT_DIR/${sample}_chr${chr}.log" 2>&1 \
      && echo "  [✓] $sample chr${chr} OK" \
      || echo "  [✗] $sample chr${chr} FALLÓ"
    ) &

    pids+=($!)

    if [[ ${#pids[@]} -ge $PARALLEL_CHRS ]]; then
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
    fi
  done

  for pid in "${pids[@]}"; do wait "$pid"; done

  # Verificar
  missing=0
  for chr in "${CHRS[@]}"; do
    [[ ! -s "$OUT_DIR/${sample}_chr${chr}_raw.vcf.gz" ]] && { echo "  [!] Falta chr${chr}"; missing=1; }
  done
  [[ "$missing" -eq 1 ]] && { echo "  [!] $sample incompleto"; continue; }

  # Merge
  echo "  [...] Mergeando $sample..."
  merged="$OUT_DIR/${sample}_merged_raw.vcf.gz"
  vcf_list=()
  for chr in "${CHRS[@]}"; do vcf_list+=("$OUT_DIR/${sample}_chr${chr}_raw.vcf.gz"); done
  bcftools concat --threads "$THREADS" -a -D -O z -o "$merged" "${vcf_list[@]}"
  bcftools index --threads "$THREADS" -t "$merged"

  # Filtrado
  echo "  [...] Filtrando..."
  filtered="$OUT_DIR/${sample}_filtered.vcf.gz"
  bcftools view --threads "$THREADS" -m2 -M2 -v snps \
    -i 'QUAL>=30 && INFO/DP>=20 && INFO/AF>=0.01 && REF!="N" && ALT!="."' \
    -O z -o "$filtered" "$merged"
  bcftools index --threads "$THREADS" -t "$filtered"
  N=$(bcftools stats "$filtered" | grep "^SN.*number of SNPs" | cut -f4)
  echo "  [✓] Filtrado: $N SNPs"

  # Liftover
  echo "  [...] Liftover..."
  unsorted="$OUT_DIR/${sample}_hg38_unsorted.vcf"
  lifted="$OUT_DIR/${sample}_hg38.vcf.gz"
  CrossMap vcf "$CHAIN" "$filtered" "$HG38" "$unsorted"
  bcftools sort --threads "$THREADS" -O z -o "$lifted" "$unsorted"
  bcftools index --threads "$THREADS" -t "$lifted"
  rm -f "$unsorted"

  # Prefijo chr
  has_chr=$(bcftools view -h "$lifted" | grep "^##contig" | head -1 | grep -c "chr" || true)
  if [[ "$has_chr" -eq 0 ]]; then
    chrmap=$(mktemp)
    for c in $(seq 1 22); do echo "$c chr$c"; done > "$chrmap"
    bcftools annotate --threads "$THREADS" --rename-chrs "$chrmap" \
      -O z -o "${lifted%.vcf.gz}_fix.vcf.gz" "$lifted"
    mv "${lifted%.vcf.gz}_fix.vcf.gz" "$lifted"
    bcftools index --threads "$THREADS" -t "$lifted"
    rm -f "$chrmap"
  fi

  N2=$(bcftools stats "$lifted" | grep "^SN.*number of SNPs" | cut -f4)
  echo "  [✓] $sample listo: $N2 SNPs en hg38"
done

echo ""
echo "✅ GATK completado. Archivos: $OUT_DIR/{SAMPLE}_hg38.vcf.gz"
