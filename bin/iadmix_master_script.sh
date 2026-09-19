#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<EOF
Uso:
  $0 --panel_name N5|N10|N11 --chr_list "21 22"|"21-22" --bam <bam> --output <prefix> [--force_build_panel] [--make_report]

Comportamiento:
- Si existe un panel precalculado en refpanels/, se usa DIRECTAMENTE (sin copiar/duplicar).
- Si no existe panel precalculado, se construye desde VCF (requiere PANEL_FILE, VCFs y bcftools).
- --make_report genera un TSV desde <prefix>.ancestry.out (requiere Rscript).
EOF
  exit 1
}

PANEL_NAME=""
CHR_LIST_RAW=""
BAM_PATH=""
OUTPUT_PREFIX=""
FORCE_BUILD_PANEL=0
MAKE_REPORT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel_name) PANEL_NAME="${2:-}"; shift 2;;
    --chr_list)   CHR_LIST_RAW="${2:-}"; shift 2;;
    --bam)        BAM_PATH="${2:-}"; shift 2;;
    --output)     OUTPUT_PREFIX="${2:-}"; shift 2;;
    --force_build_panel) FORCE_BUILD_PANEL=1; shift 1;;
    --make_report) MAKE_REPORT=1; shift 1;;
    -h|--help) usage;;
    *) echo "[X] Opción desconocida: $1" >&2; usage;;
  esac
done

[[ -n "$PANEL_NAME" && -n "$CHR_LIST_RAW" && -n "$BAM_PATH" && -n "$OUTPUT_PREFIX" ]] || usage

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export PROJ_DIR

# shellcheck disable=SC1090
source "$PROJ_DIR/env.sh"

# --- helpers ---
ts(){ date "+%F %T"; }

# log se define después; estas funciones asumen que $log ya existe
logi(){ echo "[$(ts)] [INFO] $*" | tee -a "$log"; }
loge(){ echo "[$(ts)] [X] $*" | tee -a "$log" >&2; }

expand_chr_list() {
  local s="$1"
  s="${s//,/ }"
  local out=()
  for tok in $s; do
    if [[ "$tok" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      local a="${BASH_REMATCH[1]}"; local b="${BASH_REMATCH[2]}"
      if (( a <= b )); then
        for ((c=a; c<=b; c++)); do out+=("$c"); done
      else
        for ((c=a; c>=b; c--)); do out+=("$c"); done
      fi
    else
      out+=("$tok")
    fi
  done
  printf "%s\n" "${out[@]}" | awk 'NF' | awk '!seen[$0]++'
}

# encuentra un panel precalculado con nombres variados (prioriza convenciones "FINAL_FREQUENCIES")
find_precomputed_panel() {
  local panel="$1"
  local p_low="${panel,,}"

  local candidates=(
    "$REF/FINAL_FREQUENCIES${panel}.txt"          # EJ: FINAL_FREQUENCIESN5.txt
    "$REF/FINAL_FREQUENCIES${p_low}.txt"
    "$REF/${p_low}_FINAL_FREQUENCIES.txt"         # EJ: n5_FINAL_FREQUENCIES.txt
    "$REF/${panel}_FINAL_FREQUENCIES.txt"
    "$REF/FINAL_FREQUENCIES_${panel}.txt"
    "$REF/FINAL_FREQUENCIES_${p_low}.txt"
    # Evitar usar nombres demasiado genéricos salvo que tú lo quieras
    # "$REF/${p_low}.txt"
    # "$REF/${panel}.txt"
  )

  local f
  for f in "${candidates[@]}"; do
    if [[ -s "$f" ]]; then
      echo "$f"
      return 0
    fi
  done
  return 1
}

mapfile -t CHRS < <(expand_chr_list "$CHR_LIST_RAW")
CHR_LIST="${CHRS[*]}"

POPS_FILE="$CONFIG/pops_${PANEL_NAME,,}.txt"
[[ -s "$POPS_FILE" ]] || { echo "[X] No existe POPS_FILE: $POPS_FILE" >&2; exit 1; }

# Log
log="$WORK/logs/${PANEL_NAME,,}_chr$(echo "$CHR_LIST" | tr ' ' '_').log"
mkdir -p "$(dirname "$log")"

logi "Panel: $PANEL_NAME  (POPS_FILE=$POPS_FILE)"
logi "Chr: $CHR_LIST"
logi "BAM: $BAM_PATH"
logi "OUT: $OUTPUT_PREFIX"
logi "PROJ_DIR: $PROJ_DIR"
logi "WORK: $WORK"
logi "REF: $REF"
logi "IADMIX_DIR: $IADMIX_DIR"
logi "Log: $log"

# Prechecks BAM
[[ -s "$BAM_PATH" ]] || { loge "No existe BAM: $BAM_PATH"; exit 1; }

# Acepta BAI o CSI (CSI puede aparecer en algunos entornos)
if [[ -s "${BAM_PATH}.bai" ]]; then
  logi "Índice BAM encontrado: ${BAM_PATH}.bai"
elif [[ -s "${BAM_PATH}.csi" ]]; then
  logi "Índice BAM encontrado: ${BAM_PATH}.csi"
else
  loge "Falta índice BAM (.bai o .csi): ${BAM_PATH}.bai / ${BAM_PATH}.csi"
  exit 1
fi

# Prechecks python2
command -v python2 >/dev/null 2>&1 || { loge "python2 no está disponible en PATH"; exit 1; }

# iAdmix binaries + script
for exe in ANCESTRY calculateGLL; do
  [[ -e "$IADMIX_DIR/$exe" ]] || { loge "Falta $IADMIX_DIR/$exe"; exit 1; }
  chmod +x "$IADMIX_DIR/$exe" 2>/dev/null || true
  target="$(readlink -f "$IADMIX_DIR/$exe" 2>/dev/null || true)"
  [[ -n "$target" ]] && chmod +x "$target" 2>/dev/null || true
  [[ -x "$IADMIX_DIR/$exe" ]] || { loge "No ejecutable: $IADMIX_DIR/$exe"; exit 1; }
done
[[ -s "$IADMIX_DIR/runancestry.py" ]] || { loge "Falta $IADMIX_DIR/runancestry.py"; exit 1; }

PANEL_FREQ=""

if [[ $FORCE_BUILD_PANEL -eq 0 ]]; then
  if PANEL_FREQ="$(find_precomputed_panel "$PANEL_NAME")"; then
    logi "Usando panel precalculado (sin copiar): $PANEL_FREQ"
  else
    logi "No se detectó panel precalculado para $PANEL_NAME en $REF."
  fi
fi

if [[ -z "$PANEL_FREQ" ]]; then
  # Construcción desde VCF
  [[ -s "$PANEL_FILE" ]] || { loge "Falta PANEL_FILE (requerido para construir desde VCF): $PANEL_FILE"; exit 1; }

  logi "(1/4) Listas de muestras"
  bash "$PROJ_DIR/bin/00_make_sample_lists.sh" "$POPS_FILE" | tee -a "$log"

  logi "(2/4) Frecuencias alélicas (AF TSV con bcftools)"
  bash "$PROJ_DIR/bin/01_make_af_if_needed.sh" "$POPS_FILE" "$CHR_LIST" | tee -a "$log"

  logi "(3/4) Construyendo panel FINAL_FREQUENCIES.txt"
  bash "$PROJ_DIR/bin/build_final_from_af_tsv.sh" "$POPS_FILE" "$CHR_LIST" | tee -a "$log"

  PANEL_FREQ="$REF/FINAL_FREQUENCIES.txt"
  [[ -s "$PANEL_FREQ" ]] || { loge "No se generó panel final esperado: $PANEL_FREQ"; exit 1; }
else
  logi "Saltando (1)-(3): panel ya disponible."
fi

# sanity: panel no vacío
nsnps=$(awk 'NR>1{c++} END{print c+0}' "$PANEL_FREQ")
if [[ "$nsnps" -lt 1000 ]]; then
  loge "Panel con muy pocos SNPs ($nsnps): $PANEL_FREQ"
  exit 1
fi
logi "Panel SNPs: $nsnps"

logi "(4/4) Corriendo iAdmix"
mkdir -p "$(dirname "$OUTPUT_PREFIX")"

python2 "$IADMIX_DIR/runancestry.py" \
  -f "$PANEL_FREQ" \
  --bam "$BAM_PATH" \
  -o "$OUTPUT_PREFIX" \
  --path "$IADMIX_DIR" \
  | tee -a "$log"

logi "Listo: $OUTPUT_PREFIX"

# ---------------------------
# Report (opcional) -> TSV
# ---------------------------
if [[ $MAKE_REPORT -eq 1 ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    ANCESTRY_OUT="${OUTPUT_PREFIX}.ancestry.out"
    if [[ -s "$ANCESTRY_OUT" ]]; then
      REPORT_DIR="$(dirname "$OUTPUT_PREFIX")/reports/$(basename "$OUTPUT_PREFIX")"
      mkdir -p "$REPORT_DIR"
      logi "Generando reporte TSV (R): $REPORT_DIR"

      Rscript "$PROJ_DIR/bin/03_make_report.R" \
        --input "$ANCESTRY_OUT" \
        --outdir "$REPORT_DIR" \
        --panel "$PANEL_NAME" \
        | tee -a "$log"

      logi "Reporte listo en: $REPORT_DIR"
    else
      loge "No existe: $ANCESTRY_OUT (no se pudo generar reporte)."
    fi
  else
    loge "Rscript no está disponible. Instala R o ejecuta sin --make_report."
  fi
fi

