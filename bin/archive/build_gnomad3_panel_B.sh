#!/usr/bin/env bash
# =============================================================================
# build_gnomad3_panel_B.sh
#
# Construye panel de referencia desde cero usando genotipos individuales
# del callset HGDP+TGP de gnomAD v3.1.2.
# - 100 muestras por población (balanceado)
# - AF calculadas desde genotipos reales con bcftools
# - Poblaciones: AFR, NFE, AMR, EAS, SAS
# - Panel final: refpanels/FINAL_FREQUENCIES_GNOMAD3_B.txt
#
# Uso:
#   bash bin/build_gnomad3_panel_B.sh [CHR_LIST]
#   bash bin/build_gnomad3_panel_B.sh "22"     # prueba
#   bash bin/build_gnomad3_panel_B.sh "1-22"   # completo
#
# Requisitos: bcftools
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

# ---------------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------------
POPS_ORDER=(AFR NFE AMR EAS SAS)
VCF_BASE="$PROJ_DIR/data/vcf"
SAMP_DIR="$PROJ_DIR/data/samples"
FREQ_DIR="$WORK/freq_raw_B"
COMMON_DIR="$WORK/freq_common_B"
FINAL_PANEL="$REF/FINAL_FREQUENCIES_GNOMAD3_B.txt"

mkdir -p "$FREQ_DIR" "$COMMON_DIR" "$REF"

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
echo " gnomAD v3.1.2 Panel B — 100 muestras por población"
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
# Función: calcular AF desde genotipos para 100 muestras
# ---------------------------------------------------------------------------
extract_freq_B() {
  local vcf="$1"
  local pop="$2"
  local samples_file="$3"
  local out_tsv="$4"

  if [[ -s "$out_tsv" ]]; then
    echo "[=] $(basename "$out_tsv") ya procesado."
    return
  fi

  echo "[..] Calculando AF para $pop (100 muestras)…"

  # 1) Subsetear 100 muestras
  # 2) Filtrar SNPs bialélicos PASS
  # 3) Recalcular AF desde genotipos con fill-tags
  # 4) Extraer CHROM POS REF ALT AF
  bcftools view -S "$samples_file" -f PASS -m2 -M2 -v snps "$vcf" \
    | bcftools +fill-tags -- -t AF \
    | bcftools query \
        -i 'AF > 0 && AF < 1' \
        -f "%CHROM\t%POS\t%REF\t%ALT\t%AF\n" \
    | awk 'BEGIN{OFS="\t"}
      NF==5 {
        chr=$1; sub(/^chr/,"",chr)
        pos=$2; ref=$3; alt=$4; af=$5
        if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt && af+0>0) {
          if (af < 0.0001) af = 0.0001
          if (af > 0.9999) af = 0.9999
          printf "%s\t%s\t%s\t%s\t%.6f\n", chr, pos, ref, alt, af
        }
      }' \
    > "$out_tsv"

  echo "    → $(wc -l < "$out_tsv" | tr -d ' ') SNPs"
}

# ---------------------------------------------------------------------------
# Por cromosoma
# ---------------------------------------------------------------------------
for chr in "${CHRS[@]}"; do
  echo ""
  echo "──────────────────────────────────────────"
  echo " Procesando chr${chr}…"
  echo "──────────────────────────────────────────"

  VCF_NAME="gnomad.genomes.v3.1.2.hgdp_tgp.chr${chr}.vcf.bgz"
  VCF_PATH="$VCF_BASE/$VCF_NAME"

  if [[ ! -s "$VCF_PATH" ]]; then
    echo "[!] VCF no encontrado: $VCF_PATH — descargando…"
    wget -c -q --show-progress \
      "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/${VCF_NAME}" \
      -O "$VCF_PATH"
    wget -c -q --show-progress \
      "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/${VCF_NAME}.tbi" \
      -O "${VCF_PATH}.tbi"
  fi

  # Extraer AF para cada población
  all_ok=true
  for pop in "${POPS_ORDER[@]}"; do
    SAMP="$SAMP_DIR/${pop}_samples.txt"
    OUT_TSV="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
    extract_freq_B "$VCF_PATH" "$pop" "$SAMP" "$OUT_TSV"
    if [[ ! -s "$OUT_TSV" ]]; then
      echo "[!] Sin datos para $pop chr${chr}" >&2
      all_ok=false
      break
    fi
  done

  [[ "$all_ok" == false ]] && continue

  # -------------------------------------------------------------------------
  # Intersección SNPs comunes en las 5 poblaciones
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
  # Construir filas del panel
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
echo " ✅ Panel B completado"
echo "    Archivo : $FINAL_PANEL"
echo "    SNPs    : $NSNPS"
echo "    Columnas: #chrom position rsid A1 A2 AFR NFE AMR EAS SAS"
echo "    Muestras: 100 por población (500 total)"
echo "    Genoma  : GRCh38"
echo "============================================================"
