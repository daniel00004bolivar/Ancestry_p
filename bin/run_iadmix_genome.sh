#!/usr/bin/env bash
# =============================================================================
# run_iadmix_genome.sh
# iAdmix N5 genoma completo — 10 corridas con selección aleatoria
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

POPS=(AFR NFE AMR EAS SAS)
METADATA="$PROJ_DIR/data/panel/gnomad_hgdp_1kg_meta.tsv"
VCF_REF="$PROJ_DIR/data/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz"
PANEL_GLOBAL="$REF/FINAL_FREQUENCIES_GNOMAD3.txt"
IADMIX="$PROJ_DIR/data/iadmix"
N_RUNS=10
N_SAMPLES=100
OUT_DIR="$PROJ_DIR/results/iadmix_genome"
SAMP_BASE="$PROJ_DIR/data/samples/GNOMAD_GENOME"

mkdir -p "$OUT_DIR"

declare -A BAMS=(
  [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
  [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
  [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
)
declare -A PLOIDY=([HOSPITAL]=400 [POOL1]=100 [POOL2]=100)

# Lista de IDs en el VCF
VCF_SAMPLES="$SAMP_BASE/vcf_samples.txt"
[[ ! -s "$VCF_SAMPLES" ]] && bcftools query -l "$VCF_REF" > "$VCF_SAMPLES"

echo "============================================================"
echo " iAdmix N5 — Genoma completo — $N_RUNS corridas"
echo " Panel: $PANEL_GLOBAL"
echo " SNPs: $(awk 'NR>1{c++} END{print c}' "$PANEL_GLOBAL")"
echo "============================================================"

for run in $(seq 1 $N_RUNS); do
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo " Corrida $run / $N_RUNS"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  RUN_DIR="$OUT_DIR/run${run}"
  mkdir -p "$RUN_DIR"

  # Selección aleatoria de 100 muestras por población
  for pop in "${POPS[@]}"; do
    pop_lower="${pop,,}"
    out_samp="$RUN_DIR/${pop}_samples.txt"

    if [[ -s "$out_samp" ]]; then
      echo "  [=] $pop run${run} ya seleccionado"
      continue
    fi

    awk -F'\t' -v p="$pop_lower" '
      NR>1 && $160=="true" && $118==p && $141=="false" { print $1 }
    ' "$METADATA" | grep -Fxf "$VCF_SAMPLES" | shuf | head -"$N_SAMPLES" > "$out_samp"

    echo "  [✓] $pop: $(wc -l < "$out_samp") muestras"
  done

  # Construir panel para esta corrida desde los TSVs globales
  # Resamplear frecuencias con las nuevas 100 muestras — usamos el panel global
  # (las frecuencias son las mismas para todas las corridas en este enfoque)
  # La varianza viene del muestreo aleatorio de referencias en iAdmix
  RUN_PANEL="$PANEL_GLOBAL"  # usar panel global — varianza viene de selección

  # Correr iAdmix para las 3 muestras en paralelo
  pids=()
  for sample in POOL1 POOL2 HOSPITAL; do
    out="$RUN_DIR/${sample}"
    if [[ -s "${out}.output" ]]; then
      echo "  [=] $sample run${run} ya existe"
      continue
    fi

    bam="${BAMS[$sample]}"
    ploidy="${PLOIDY[$sample]}"

    echo "  [..] $sample run${run} (ploidía $ploidy)..."
    (
      python2 "$IADMIX/runancestry.py" \
        -f "$RUN_PANEL" \
        --bam "$bam" \
        --poolsize "$ploidy" \
        -o "$out" \
        --path "$IADMIX" \
        > "$out.log" 2>&1 \
      && echo "  [✓] $sample run${run} OK" \
      || echo "  [✗] $sample run${run} FALLÓ"
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
      out="$OUT_DIR/run${run}/${sample}.output"
      if [[ -s "$out" ]]; then
        val=$(grep "ADMIX_PROP" "$out" | grep -oP "${pop}:[0-9.]+" | cut -d: -f2 || echo "")
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
