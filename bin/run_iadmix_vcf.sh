#!/usr/bin/env bash
# =============================================================================
# run_iadmix_vcf.sh
# iAdmix N5 genoma completo — usa VCFs hg38 en lugar de BAMs
# 10 corridas con el mismo panel (varianza viene de selección de referencias)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

POPS=(AFR NFE AMR EAS SAS)
PANEL_GLOBAL="$REF/FINAL_FREQUENCIES_GNOMAD3.txt"
IADMIX="$PROJ_DIR/data/iadmix"
N_RUNS=10
OUT_DIR="$PROJ_DIR/results/iadmix_genome"
GATK_DIR="$WORK/gatk_genome"

mkdir -p "$OUT_DIR"

declare -A VCFS=(
  [HOSPITAL]="$GATK_DIR/HOSPITAL_hg38.vcf.gz"
  [POOL1]="$GATK_DIR/POOL1_hg38.vcf.gz"
  [POOL2]="$GATK_DIR/POOL2_hg38.vcf.gz"
)
declare -A POOLSIZE=([HOSPITAL]=400 [POOL1]=100 [POOL2]=100)

echo "============================================================"
echo " iAdmix N5 — Genoma completo — $N_RUNS corridas"
echo " Input: VCFs hg38"
echo " Panel: $PANEL_GLOBAL"
echo " SNPs panel: $(awk 'NR>1{c++} END{print c}' "$PANEL_GLOBAL")"
echo "============================================================"

# Verificar que los VCFs existen
for sample in "${!VCFS[@]}"; do
  vcf="${VCFS[$sample]}"
  [[ -s "$vcf" ]] || { echo "[X] VCF no encontrado: $vcf" >&2; exit 1; }
  echo "  ✓ $sample: $(basename $vcf)"
done

for run in $(seq 1 $N_RUNS); do
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo " Corrida $run / $N_RUNS"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  RUN_DIR="$OUT_DIR/run${run}"
  mkdir -p "$RUN_DIR"

  # Correr iAdmix para las 3 muestras en paralelo
  pids=()
  for sample in POOL1 POOL2 HOSPITAL; do
    out="$RUN_DIR/${sample}"
    if [[ -s "${out}.ancestry.out" ]]; then
      echo "  [=] $sample run${run} ya existe"
      # Mostrar resultado
      grep "ADMIX_PROP" "${out}.ancestry.out" | tail -1 | \
        grep -oP "(AFR|NFE|AMR|EAS|SAS):[0-9.]+" | tr '\n' ' '
      echo ""
      continue
    fi

    vcf="${VCFS[$sample]}"
    poolsize="${POOLSIZE[$sample]}"

    echo "  [..] $sample run${run} (poolsize $poolsize)..."
    (
      python2 "$IADMIX/runancestry.py" \
        -f "$PANEL_GLOBAL" \
        --vcf "$vcf" \
        -p "$poolsize" \
        -o "$out" \
        --path "$IADMIX" \
        > "$out.log" 2>&1 \
      && echo "  [✓] $sample run${run}: $(grep 'ADMIX_PROP' ${out}.ancestry.out | tail -1 | grep -oP '(AFR|NFE|AMR|EAS|SAS):[0-9.]+' | tr '\n' ' ')" \
      || echo "  [✗] $sample run${run} FALLÓ — ver $out.log"
    ) &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid"; done

done

# Consolidar resultados
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " Consolidando resultados..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

SUMMARY="$OUT_DIR/summary_N5_genome.tsv"
printf "Sample\tPop\tMean\tSD\tMin\tMax\n" > "$SUMMARY"

for sample in POOL1 POOL2 HOSPITAL; do
  for pop in "${POPS[@]}"; do
    vals=()
    for run in $(seq 1 $N_RUNS); do
      out="$OUT_DIR/run${run}/${sample}.ancestry.out"
      if [[ -s "$out" ]]; then
        val=$(grep "ADMIX_PROP" "$out" | tail -1 | \
              grep -oP "${pop}:[0-9.]+" | cut -d: -f2 || echo "")
        [[ -n "$val" ]] && vals+=("$val")
      fi
    done

    if [[ ${#vals[@]} -gt 0 ]]; then
      printf '%s\n' "${vals[@]}" | awk -v s="$sample" -v p="$pop" '
        BEGIN{n=0;sum=0;sum2=0;min=999;max=-999}
        {n++;sum+=$1;sum2+=$1*$1;if($1<min)min=$1;if($1>max)max=$1}
        END{
          mean=sum/n
          sd=sqrt((sum2-n*mean*mean)/(n>1?n-1:1))
          printf "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",s,p,mean,sd,min,max
        }' >> "$SUMMARY"
    fi
  done
done

echo ""
echo "============================================================"
echo " ✅ iAdmix completado"
echo " Resumen: $SUMMARY"
echo "============================================================"
cat "$SUMMARY"
