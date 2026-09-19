#!/usr/bin/env bash
# =============================================================================
# run_bootstrap.sh
#
# Corre el pipeline completo N veces con semillas diferentes para estimar
# la sensibilidad de los resultados a la selección de muestras de referencia.
#
# Para cada iteración:
#   1. Construye panel con semilla distinta (muestras aleatorias diferentes)
#   2. Corre IAdmix para cada BAM
#   3. Recopila resultados
# Al final calcula media ± SD por población y por muestra.
#
# Uso:
#   bash bin/run_bootstrap.sh [OPCIONES]
#
# Opciones:
#   --bams    "HOSPITAL,POOL1,POOL2"   nombres de las muestras (default: todas)
#   --chrs    "22"                     cromosomas (default: 22)
#   --pops    "AFR,NFE,AMR,EAS,SAS"   poblaciones (default: todas)
#   --n       100                      muestras por población (default: 100)
#   --threads 20                       cores (default: 4)
#   --iters   10                       número de iteraciones (default: 10)
#   --out     "BOOTSTRAP"              nombre del experimento (default: BOOTSTRAP)
#
# Ejemplo:
#   bash bin/run_bootstrap.sh --chrs "22" --pops "AFR,NFE,AMR,EAS,SAS" --threads 20 --iters 10
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
CHR_INPUT="22"
POPS_INPUT="AFR,NFE,AMR,EAS,SAS"
N_SAMPLES=100
THREADS=4
N_ITERS=10
EXP_NAME="BOOTSTRAP"

# BAMs disponibles — nombre:ruta
declare -A BAM_MAP=(
  [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
  [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
  [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
)
BAM_NAMES_DEFAULT="HOSPITAL,POOL1,POOL2"
BAM_NAMES="$BAM_NAMES_DEFAULT"

# ---------------------------------------------------------------------------
# Parse argumentos
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bams)    BAM_NAMES="$2";   shift 2 ;;
    --chrs)    CHR_INPUT="$2";   shift 2 ;;
    --pops)    POPS_INPUT="$2";  shift 2 ;;
    --n)       N_SAMPLES="$2";   shift 2 ;;
    --threads) THREADS="$2";     shift 2 ;;
    --iters)   N_ITERS="$2";     shift 2 ;;
    --out)     EXP_NAME="$2";    shift 2 ;;
    *) echo "[X] Argumento desconocido: $1" >&2; exit 1 ;;
  esac
done

IFS=',' read -ra POPS_ORDER <<< "$POPS_INPUT"
IFS=',' read -ra BAM_LIST   <<< "$BAM_NAMES"

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
BOOT_DIR="$PROJ_DIR/work/bootstrap_${EXP_NAME}"
RESULTS_DIR="$BOOT_DIR/results"
mkdir -p "$RESULTS_DIR"

# Archivo con todos los resultados crudos
RAW_RESULTS="$BOOT_DIR/raw_results.tsv"
SUMMARY="$BOOT_DIR/summary_mean_sd.tsv"

echo "============================================================"
echo " Bootstrap — ${EXP_NAME}"
echo " Iteraciones : ${N_ITERS}"
echo " Cromosomas  : ${CHR_INPUT}"
echo " Poblaciones : ${POPS_ORDER[*]}"
echo " Muestras    : ${N_SAMPLES} por población"
echo " BAMs        : ${BAM_LIST[*]}"
echo " Cores       : ${THREADS}"
echo " Resultados  : ${BOOT_DIR}"
echo "============================================================"

# Header del archivo de resultados crudos
{
  printf "iter\tseed\tbam"
  for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done
  printf "\n"
} > "$RAW_RESULTS"

# ---------------------------------------------------------------------------
# Loop de iteraciones
# ---------------------------------------------------------------------------
for iter in $(seq 1 "$N_ITERS"); do
  SEED=$iter
  PANEL_NAME="${EXP_NAME}_seed${SEED}"

  echo ""
  echo "════════════════════════════════════════════"
  echo " Iteración ${iter}/${N_ITERS} — semilla ${SEED}"
  echo "════════════════════════════════════════════"

  # Construir panel con esta semilla
  bash "$SCRIPT_DIR/build_gnomad3_panel.sh" \
    --chrs  "$CHR_INPUT" \
    --pops  "$POPS_INPUT" \
    --n     "$N_SAMPLES" \
    --threads "$THREADS" \
    --seed  "$SEED" \
    --out   "$PANEL_NAME"

  # Crear config de poblaciones para IAdmix
  POPS_CFG="$PROJ_DIR/config/pops_${PANEL_NAME,,}.txt"
  printf "%s\n" "${POPS_ORDER[@]}" > "$POPS_CFG"

  # Correr IAdmix para cada BAM
  for bam_name in "${BAM_LIST[@]}"; do
    BAM_PATH="${BAM_MAP[$bam_name]}"
    OUTPUT="$RESULTS_DIR/${bam_name}_seed${SEED}.output"

    echo "  [..] IAdmix: ${bam_name} (semilla ${SEED})…"

    bash "$SCRIPT_DIR/iadmix_master_script.sh" \
      --panel_name "$PANEL_NAME" \
      --chr_list   "$CHR_INPUT" \
      --bam        "$BAM_PATH" \
      --output     "$OUTPUT"

    # Extraer proporciones del output
    PROPS_LINE=$(grep "^final maxval" "${OUTPUT}.ancestry.out" | tail -1)

    # Parsear proporciones: AFR:0.76 NFE:0.00 ...
    ROW="${iter}\t${SEED}\t${bam_name}"
    for pop in "${POPS_ORDER[@]}"; do
      VAL=$(echo "$PROPS_LINE" | grep -oP "${pop}:\K[0-9.]+" || echo "NA")
      ROW="${ROW}\t${VAL}"
    done
    printf "%b\n" "$ROW" >> "$RAW_RESULTS"

    echo "  [✓] ${bam_name} semilla ${SEED} completado"
  done

  echo "  [✓] Iteración ${iter} completada"
done

# ---------------------------------------------------------------------------
# Calcular media y SD por BAM y por población
# ---------------------------------------------------------------------------
echo ""
echo "Calculando estadísticas…"

{
  printf "bam\tstat"
  for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done
  printf "\n"

  for bam_name in "${BAM_LIST[@]}"; do
    # Extraer columnas de proporciones para este BAM
    awk -F'\t' -v bam="$bam_name" -v pops="${POPS_ORDER[*]}" -v np="${#POPS_ORDER[@]}" '
      BEGIN {
        split(pops, P, " ")
        for (i=1; i<=np; i++) { sum[i]=0; sum2[i]=0; n[i]=0 }
      }
      NR>1 && $3==bam {
        for (i=1; i<=np; i++) {
          val = $(3+i)
          if (val != "NA") {
            sum[i]  += val
            sum2[i] += val*val
            n[i]++
          }
        }
      }
      END {
        # Media
        printf "%s\tmean", bam
        for (i=1; i<=np; i++) {
          if (n[i]>0) printf "\t%.4f", sum[i]/n[i]
          else printf "\tNA"
        }
        printf "\n"

        # SD
        printf "%s\tsd", bam
        for (i=1; i<=np; i++) {
          if (n[i]>1) {
            mean = sum[i]/n[i]
            var  = (sum2[i] - n[i]*mean*mean) / (n[i]-1)
            printf "\t%.4f", (var>0 ? sqrt(var) : 0)
          } else printf "\tNA"
        }
        printf "\n"
      }
    ' "$RAW_RESULTS"
  done
} > "$SUMMARY"

echo ""
echo "============================================================"
echo " ✅ Bootstrap completado"
echo "    Resultados crudos : $RAW_RESULTS"
echo "    Resumen media±SD  : $SUMMARY"
echo "============================================================"
echo ""
echo "Resumen:"
cat "$SUMMARY"
