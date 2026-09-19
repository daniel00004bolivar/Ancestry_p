#!/usr/bin/env bash
# =============================================================================
# build_gnomad3_panel_v2.sh
#
# Construye un panel de referencia desde cero usando gnomAD v3.1.2 sites VCF.
# - Sin liftover
# - Sin dependencia del panel N5
# - AF calculadas directamente desde AC_pop/AN_pop del INFO field
# - Poblaciones: AFR, NFE, AMR, EAS, SAS
# - Panel final: refpanels/FINAL_FREQUENCIES_GNOMAD3.txt
#
# Uso:
#   bash bin/build_gnomad3_panel_v2.sh [CHR_LIST]
#   bash bin/build_gnomad3_panel_v2.sh "22"       # prueba
#   bash bin/build_gnomad3_panel_v2.sh "1-22"     # completo
#
# Requisitos: bcftools, wget
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

# ---------------------------------------------------------------------------
# Configuración de poblaciones
# NFE = Non-Finnish European (nombre correcto en gnomAD, no EUR)
# ---------------------------------------------------------------------------
declare -A POP_AC=(
  [AFR]="AC_afr"
  [NFE]="AC_nfe"
  [AMR]="AC_amr"
  [EAS]="AC_eas"
  [SAS]="AC_sas"
)
declare -A POP_AN=(
  [AFR]="AN_afr"
  [NFE]="AN_nfe"
  [AMR]="AN_amr"
  [EAS]="AN_eas"
  [SAS]="AN_sas"
)
POPS_ORDER=(AFR NFE AMR EAS SAS)

GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"
FINAL_PANEL="$REF/FINAL_FREQUENCIES_GNOMAD3.txt"
FREQ_DIR="$WORK/freq_raw"
COMMON_DIR="$WORK/freq_common"
VCF_DIR="$PROJ_DIR/data/vcf"

mkdir -p "$FREQ_DIR" "$COMMON_DIR" "$VCF_DIR" "$REF"

# ---------------------------------------------------------------------------
# Parse CHR_LIST
# ---------------------------------------------------------------------------
RAW_CHR="${1:-1-22}"
CHRS=()
for token in $RAW_CHR; do
  if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    for ((c=${BASH_REMATCH[1]}; c<=${BASH_REMATCH[2]}; c++)); do CHRS+=("$c"); done
  else
    CHRS+=("$token")
  fi
done

echo "============================================================"
echo " gnomAD v3.1.2 panel builder — desde cero"
echo " Cromosomas : ${CHRS[*]}"
echo " Poblaciones: ${POPS_ORDER[*]}"
echo " Panel final: $FINAL_PANEL"
echo "============================================================"

# ---------------------------------------------------------------------------
# Header del panel final
# ---------------------------------------------------------------------------
{
  printf "#chrom\tposition\trsid\tA1\tA2"
  for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done
  printf "\n"
} > "$FINAL_PANEL"

# ---------------------------------------------------------------------------
# Función: descargar VCF + índice si no existe
# ---------------------------------------------------------------------------
download_vcf() {
  local vcf_path="$1"
  local vcf_url="$2"
  if [[ ! -s "$vcf_path" ]]; then
    echo "[↓] Descargando $(basename "$vcf_path")…"
    wget -c -q --show-progress -O "$vcf_path" "$vcf_url"
  else
    echo "[=] $(basename "$vcf_path") ya en caché."
  fi
  if [[ ! -s "${vcf_path}.tbi" ]]; then
    echo "[↓] Descargando índice .tbi…"
    wget -c -q --show-progress -O "${vcf_path}.tbi" "${vcf_url}.tbi"
  fi
}

# ---------------------------------------------------------------------------
# Función: extraer AF desde sites VCF usando AC_pop/AN_pop
# Una sola pasada lineal por el VCF — sin -R, sin liftover
# Output: CHROM POS REF ALT AF  (sin prefijo chr, posiciones GRCh38)
# ---------------------------------------------------------------------------
extract_freq() {
  local vcf="$1"
  local ac_tag="$2"
  local an_tag="$3"
  local out_tsv="$4"

  if [[ -s "$out_tsv" ]]; then
    echo "[=] $(basename "$out_tsv") ya procesado."
    return
  fi

  echo "[..] Extrayendo ${ac_tag}/${an_tag}…"

  bcftools view -f PASS -m2 -M2 -v snps "$vcf" \
    | bcftools query \
        -i "${ac_tag} != \".\" && ${an_tag} > 0" \
        -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/${ac_tag}\t%INFO/${an_tag}\n" \
    | awk 'BEGIN{OFS="\t"}
      NF==6 {
        chr=$1; sub(/^chr/,"",chr)
        pos=$2; ref=$3; alt=$4; ac=$5; an=$6
        if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt && an+0>0 && ac+0>=0) {
          af = ac/an
          if (af < 0.0001) af = 0.0001
          if (af > 0.9999) af = 0.9999
          printf "%s\t%s\t%s\t%s\t%.6f\n", chr, pos, ref, alt, af
        }
      }' \
    > "$out_tsv"

  echo "    → $(wc -l < "$out_tsv" | tr -d ' ') SNPs"
}

# ---------------------------------------------------------------------------
# Por cromosoma: descargar, extraer, construir panel
# ---------------------------------------------------------------------------
for chr in "${CHRS[@]}"; do
  echo ""
  echo "──────────────────────────────────────────"
  echo " Procesando chr${chr}…"
  echo "──────────────────────────────────────────"

  VCF_NAME="gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz"
  VCF_PATH="$VCF_DIR/$VCF_NAME"
  VCF_URL="${GNOMAD_BASE}/${VCF_NAME}"

  download_vcf "$VCF_PATH" "$VCF_URL"

  # Extraer AF para cada población
  all_ok=true
  for pop in "${POPS_ORDER[@]}"; do
    OUT_TSV="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
    extract_freq "$VCF_PATH" "${POP_AC[$pop]}" "${POP_AN[$pop]}" "$OUT_TSV"
    if [[ ! -s "$OUT_TSV" ]]; then
      echo "[!] Sin datos para $pop chr${chr}" >&2
      all_ok=false
      break
    fi
  done

  # Borrar VCF para liberar espacio
  echo "[🗑] Liberando VCF chr${chr}…"
  rm -f "$VCF_PATH" "${VCF_PATH}.tbi"

  [[ "$all_ok" == false ]] && continue

  # -------------------------------------------------------------------------
  # Intersección de SNPs presentes en las 5 poblaciones
  # -------------------------------------------------------------------------
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

  N_POPS=${#POPS_ORDER[@]}
  COMMON_IDS="$COMMON_DIR/chr${chr}.ids"

  cat "${id_files[@]}" | sort | awk -v n="$N_POPS" '
    { if ($0==prev){c++} else { if (prev!="" && c==n) print prev; prev=$0; c=1 } }
    END{ if (prev!="" && c==n) print prev }
  ' > "$COMMON_IDS"

  N_COMMON=$(wc -l < "$COMMON_IDS" | tr -d ' ')
  echo "[✓] SNPs comunes chr${chr}: ${N_COMMON}"
  [[ "$N_COMMON" -eq 0 ]] && continue

  # -------------------------------------------------------------------------
  # Construir filas del panel final
  # -------------------------------------------------------------------------
  files=("$COMMON_IDS")
  for pop in "${POPS_ORDER[@]}"; do files+=("$COMMON_DIR/${pop}_chr${chr}.tab"); done

  awk -v pops_str="$(IFS=' '; echo "${POPS_ORDER[*]}")" '
    BEGIN{
      OFS="\t"
      split(pops_str, P, " ")
      np=length(P)
    }
    ARGIND==1 { ids[$1]=1; next }
    ARGIND>1  {
      key=$1; chr=$2; pos=$3; ref=$4; alt=$5; af=$6
      idx=ARGIND-1; pop=P[idx]
      freq[pop,key]=af
      if(!(key in meta)) {
        rsid="rs"chr"_"pos"_"ref"_"alt
        meta[key]=chr OFS pos OFS rsid OFS ref OFS alt
      }
      next
    }
    END{
      for(k in ids){
        if(!(k in meta)) continue
        printf "%s", meta[k]
        for(i=1; i<=np; i++){
          pop=P[i]
          val=((pop SUBSEP k) in freq) ? freq[pop,k] : "0.000100"
          printf OFS "%.6f", val+0
        }
        printf "\n"
      }
    }' "${files[@]}" \
    | sort -k1,1n -k2,2n \
    >> "$FINAL_PANEL"

  echo "[✓] chr${chr} añadido al panel."
done

# ---------------------------------------------------------------------------
# Validación final
# ---------------------------------------------------------------------------
NF_EXPECTED=$((5 + ${#POPS_ORDER[@]}))
awk -v nf="$NF_EXPECTED" 'BEGIN{FS=OFS="\t"} NR==1{print;next} NF==nf{print}' \
  "$FINAL_PANEL" > "${FINAL_PANEL}.tmp" && mv "${FINAL_PANEL}.tmp" "$FINAL_PANEL"

NSNPS=$(awk 'NR>1{c++} END{print c+0}' "$FINAL_PANEL")

echo ""
echo "============================================================"
echo " ✅ Panel gnomAD v3.1.2 completado"
echo "    Archivo : $FINAL_PANEL"
echo "    SNPs    : $NSNPS"
echo "    Columnas: #chrom position rsid A1 A2 AFR NFE AMR EAS SAS"
echo "    Genoma  : GRCh38"
echo "============================================================"
