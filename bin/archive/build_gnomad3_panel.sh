#!/usr/bin/env bash
# =============================================================================
# build_gnomad3_panel.sh  — v4.1
#
# Uso:
#   bash bin/build_gnomad3_panel.sh [OPCIONES]
#
# Opciones:
#   --chrs    "22"                    cromosomas (default: 1-22)
#   --pops    "AFR,NFE,AMR,EAS,SAS"  poblaciones (default: todas)
#   --n       100                     muestras por población (default: 100)
#   --threads 20                      cores totales (default: 4)
#   --min_af  0.01                    filtro mínimo AF (default: 0.01)
#   --out     "GNOMAD3"               nombre del panel (default: GNOMAD3)
#   --seed    42                      semilla aleatoria (default: 42)
#
# Ejemplos:
#   bash bin/build_gnomad3_panel.sh --chrs "22" --pops "AFR,NFE,AMR,EAS,SAS" --threads 20
#   bash bin/build_gnomad3_panel.sh --chrs "1-22" --pops "AFR,NFE,AMR,EAS,SAS" --threads 20 --n 150 --seed 99
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
CHR_INPUT="1-22"
POPS_INPUT="AFR,NFE,AMR,EAS,SAS"
N_SAMPLES=100
THREADS=4
MIN_AF=0.01
PANEL_NAME="GNOMAD3"
SEED=42

# ---------------------------------------------------------------------------
# Parse argumentos
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chrs)    CHR_INPUT="$2";    shift 2 ;;
    --pops)    POPS_INPUT="$2";   shift 2 ;;
    --n)       N_SAMPLES="$2";    shift 2 ;;
    --threads) THREADS="$2";      shift 2 ;;
    --min_af)  MIN_AF="$2";       shift 2 ;;
    --out)     PANEL_NAME="$2";   shift 2 ;;
    --seed)    SEED="$2";         shift 2 ;;
    *) echo "[X] Argumento desconocido: $1" >&2; exit 1 ;;
  esac
done

IFS=',' read -ra POPS_ORDER <<< "$POPS_INPUT"
N_POPS=${#POPS_ORDER[@]}
THREADS_PER_POP=$(( THREADS / N_POPS ))
[[ "$THREADS_PER_POP" -lt 1 ]] && THREADS_PER_POP=1

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
VCF_BASE="$PROJ_DIR/data/vcf"
SAMP_DIR="$PROJ_DIR/data/samples/${PANEL_NAME}"
FREQ_DIR="$WORK/freq_raw_${PANEL_NAME}"
COMMON_DIR="$WORK/freq_common_${PANEL_NAME}"
FINAL_PANEL="$REF/FINAL_FREQUENCIES_${PANEL_NAME}.txt"
METADATA="$PROJ_DIR/data/panel/gnomad_hgdp_1kg_meta.tsv"

mkdir -p "$SAMP_DIR" "$FREQ_DIR" "$COMMON_DIR" "$REF"

# ---------------------------------------------------------------------------
# Parse cromosomas
# ---------------------------------------------------------------------------
CHRS=()
for token in $CHR_INPUT; do
  if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    for ((c=${BASH_REMATCH[1]}; c<=${BASH_REMATCH[2]}; c++)); do CHRS+=("$c"); done
  else
    CHRS+=("$token")
  fi
done

echo "============================================================"
echo " gnomAD v3.1.2 Panel — ${PANEL_NAME}"
echo " Cromosomas  : ${CHRS[*]}"
echo " Poblaciones : ${POPS_ORDER[*]}"
echo " Muestras    : ${N_SAMPLES} por población"
echo " Cores       : ${THREADS} (${THREADS_PER_POP} por pop)"
echo " Filtro AF   : >= ${MIN_AF}"
echo " Semilla     : ${SEED}"
echo " Panel final : ${FINAL_PANEL}"
echo "============================================================"

# ---------------------------------------------------------------------------
# PASO 1: Selección de muestras
# Columnas del metadata:
#   $1   = sample_id
#   $118 = population_inference.pop  (afr, nfe, amr, eas, sas)
#   $160 = high_quality              (true/false)
#   $141 = sample_filters.release_related (true/false)
# ---------------------------------------------------------------------------
echo ""
echo "[1/4] Seleccionando muestras…"

VCF_REF="$VCF_BASE/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz"
[[ ! -s "$VCF_REF" ]] && { echo "[X] VCF chr22 no encontrado: $VCF_REF" >&2; exit 1; }

VCF_SAMPLES="$SAMP_DIR/vcf_samples.txt"
if [[ ! -s "$VCF_SAMPLES" ]]; then
  echo "    Extrayendo IDs del VCF…"
  bcftools query -l "$VCF_REF" > "$VCF_SAMPLES"
fi
echo "    Muestras en VCF: $(wc -l < "$VCF_SAMPLES")"

for pop in "${POPS_ORDER[@]}"; do
  pop_lower="${pop,,}"
  OUT_SAMP="$SAMP_DIR/${pop}_samples.txt"

  if [[ -s "$OUT_SAMP" ]]; then
    echo "    [=] $pop: ya seleccionado ($(wc -l < "$OUT_SAMP") muestras)"
    continue
  fi

  echo "    Seleccionando $pop…"

  # Filtrar metadata:
  # - población correcta ($118)
  # - high_quality=true ($160)
  # - NO relacionado: release_related=false ($141)
  # - presente en el VCF (hash en awk)
  awk -F'\t' -v p="$pop_lower" -v vcf="$VCF_SAMPLES" '
    BEGIN { while ((getline line < vcf) > 0) vcf_ids[line] = 1 }
    NR > 1 &&
    $118 == p &&
    $160 == "true" &&
    $141 == "false" &&
    ($1 in vcf_ids) { print $1 }
  ' "$METADATA" | shuf --random-source=<(yes "$SEED") | head -"$N_SAMPLES" > "$OUT_SAMP"

  N_SEL=$(wc -l < "$OUT_SAMP")

  # Si no hay suficientes, relajar filtro de relacionados
  if [[ "$N_SEL" -eq 0 ]]; then
    echo "    [!] Sin muestras para $pop — relajando filtro de relacionados" >&2
    awk -F'\t' -v p="$pop_lower" -v vcf="$VCF_SAMPLES" '
      BEGIN { while ((getline line < vcf) > 0) vcf_ids[line] = 1 }
      NR > 1 && $118 == p && $160 == "true" && ($1 in vcf_ids) { print $1 }
    ' "$METADATA" | shuf --random-source=<(yes "$SEED") | head -"$N_SAMPLES" > "$OUT_SAMP"
    N_SEL=$(wc -l < "$OUT_SAMP")
  fi

  echo "    [✓] $pop: $N_SEL muestras seleccionadas"
done

# ---------------------------------------------------------------------------
# PASO 2: Calcular AF por población en paralelo
# ---------------------------------------------------------------------------
echo ""
echo "[2/4] Calculando frecuencias alélicas…"

extract_freq() {
  local vcf="$1"
  local pop="$2"
  local samples_file="$3"
  local out_tsv="$4"
  local min_af="$5"
  local threads="$6"
  local chr="$7"

  if [[ -s "$out_tsv" ]]; then
    echo "    [=] $pop chr${chr} ya procesado."
    return
  fi

  local max_af
  max_af=$(awk "BEGIN{printf \"%.4f\", 1 - $min_af}")

  echo "    [..] $pop chr${chr} iniciando…"

  bcftools view --threads "$threads" -S "$samples_file" \
      -f PASS -m2 -M2 -v snps "$vcf" \
    | bcftools +fill-tags -- -t AF \
    | bcftools query \
        -i "AF >= ${min_af} && AF <= ${max_af}" \
        -f "%CHROM\t%POS\t%REF\t%ALT\t%AF\n" \
    | awk -v maf="$min_af" 'BEGIN{OFS="\t"}
      NF==5 {
        chr=$1; sub(/^chr/,"",chr)
        pos=$2; ref=$3; alt=$4; af=$5
        if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt && af+0>=maf+0) {
          if (af < 0.0001) af = 0.0001
          if (af > 0.9999) af = 0.9999
          printf "%s\t%s\t%s\t%s\t%.6f\n", chr, pos, ref, alt, af
        }
      }' > "$out_tsv"

  echo "    [✓] $pop chr${chr}: $(wc -l < "$out_tsv" | tr -d ' ') SNPs"
}

export -f extract_freq

# ---------------------------------------------------------------------------
# PASO 3 y 4: Por cromosoma — intersección + panel final
# ---------------------------------------------------------------------------
{
  printf "#chrom\tposition\trsid\tA1\tA2"
  for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done
  printf "\n"
} > "$FINAL_PANEL"

for chr in "${CHRS[@]}"; do
  echo ""
  echo "──────────────────────────────────────────"
  echo " chr${chr}"
  echo "──────────────────────────────────────────"

  VCF_NAME="gnomad.genomes.v3.1.2.hgdp_tgp.chr${chr}.vcf.bgz"
  VCF_PATH="$VCF_BASE/$VCF_NAME"

  if [[ ! -s "$VCF_PATH" ]]; then
    echo "[↓] Descargando $VCF_NAME…"
    wget -c -q --show-progress \
      "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/${VCF_NAME}" \
      -O "$VCF_PATH"
    wget -c -q --show-progress \
      "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/${VCF_NAME}.tbi" \
      -O "${VCF_PATH}.tbi"
  fi

  # Lanzar todas las poblaciones en paralelo
  pids=()
  for pop in "${POPS_ORDER[@]}"; do
    SAMP="$SAMP_DIR/${pop}_samples.txt"
    OUT_TSV="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
    extract_freq "$VCF_PATH" "$pop" "$SAMP" "$OUT_TSV" \
      "$MIN_AF" "$THREADS_PER_POP" "$chr" &
    pids+=($!)
  done

  all_ok=true
  for pid in "${pids[@]}"; do
    wait "$pid" || { all_ok=false; }
  done

  for pop in "${POPS_ORDER[@]}"; do
    if [[ ! -s "$FREQ_DIR/${pop}_chr${chr}_af.tsv" ]]; then
      echo "[!] Sin datos para $pop chr${chr}" >&2
      all_ok=false
    fi
  done
  [[ "$all_ok" == false ]] && continue

  # Intersección SNPs comunes en todas las poblaciones
  id_files=()
  for pop in "${POPS_ORDER[@]}"; do
    IN="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
    TAB="$COMMON_DIR/${pop}_chr${chr}.tab"
    IDS="$COMMON_DIR/${pop}_chr${chr}.ids"

    awk 'BEGIN{OFS="\t"}
      NF>=5 {
        chr=$1; pos=$2; ref=$3; alt=$4; af=$5
        if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt) {
          key=chr ":" pos ":" ref ":" alt
          print key, chr, pos, ref, alt, af
        }
      }' "$IN" > "$TAB"

    cut -f1 "$TAB" | sort -u > "$IDS"
    id_files+=("$IDS")
  done

  COMMON_IDS="$COMMON_DIR/chr${chr}.ids"
  cat "${id_files[@]}" | sort | awk -v n="$N_POPS" '
    { if ($0==prev){c++} else { if (prev!="" && c==n) print prev; prev=$0; c=1 } }
    END{ if (prev!="" && c==n) print prev }
  ' > "$COMMON_IDS"

  N_COMMON=$(wc -l < "$COMMON_IDS" | tr -d ' ')
  echo "[3/4] SNPs comunes chr${chr}: ${N_COMMON}"
  [[ "$N_COMMON" -eq 0 ]] && continue

  # Construir filas del panel
  files=("$COMMON_IDS")
  for pop in "${POPS_ORDER[@]}"; do files+=("$COMMON_DIR/${pop}_chr${chr}.tab"); done

  awk -v pops_str="$(IFS=' '; echo "${POPS_ORDER[*]}")" '
    BEGIN{ OFS="\t"; split(pops_str,P," "); np=length(P) }
    ARGIND==1 { ids[$1]=1; next }
    ARGIND>1  {
      key=$1; chr=$2; pos=$3; ref=$4; alt=$5; af=$6
      idx=ARGIND-1; pop=P[idx]
      freq[pop,key]=af
      if(!(key in meta)){
        meta[key]=chr OFS pos OFS "rs"chr"_"pos"_"ref"_"alt OFS ref OFS alt
      }
      next
    }
    END{
      for(k in ids){
        if(!(k in meta)) continue
        printf "%s", meta[k]
        for(i=1;i<=np;i++){
          pop=P[i]
          val=((pop SUBSEP k) in freq) ? freq[pop,k] : 0.010000
          printf OFS "%.6f", val+0
        }
        printf "\n"
      }
    }' "${files[@]}" | sort -k1,1n -k2,2n >> "$FINAL_PANEL"

  echo "[4/4] chr${chr} añadido al panel."
done

# ---------------------------------------------------------------------------
# Validación final
# ---------------------------------------------------------------------------
NF_EXPECTED=$((5 + N_POPS))
awk -v nf="$NF_EXPECTED" 'BEGIN{FS=OFS="\t"} NR==1{print;next} NF==nf{print}' \
  "$FINAL_PANEL" > "${FINAL_PANEL}.tmp" && mv "${FINAL_PANEL}.tmp" "$FINAL_PANEL"

NSNPS=$(awk 'NR>1{c++} END{print c+0}' "$FINAL_PANEL")

echo ""
echo "============================================================"
echo " ✅ Panel completado: ${PANEL_NAME}"
echo "    Archivo    : $FINAL_PANEL"
echo "    SNPs       : $NSNPS"
echo "    Poblaciones: ${POPS_ORDER[*]}"
echo "    Muestras   : ${N_SAMPLES} por población (aleatorio, semilla ${SEED}, sin relacionados)"
echo "    Filtro AF  : >= ${MIN_AF}"
echo "    Genoma     : GRCh38"
echo "============================================================"
