#!/usr/bin/env bash
# =============================================================================
# 02_build_panel_genome.sh
#
# Construye panel N5 (AFR, NFE, AMR, EAS, SAS) para genoma completo.
#
# ESTRATEGIA DE DISCO:
#   Solo 140 GB libres. Cada VCF de gnomAD pesa ~30-60 GB.
#   → Descargamos un cromosoma, extraemos frecuencias, borramos el VCF.
#   → Acumulamos el panel en un archivo único.
#
# Uso: bash 02_build_panel_genome.sh [THREADS]
# =============================================================================
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$PROJ/env.sh"

THREADS="${1:-20}"
THREADS_PER_POP=$(( THREADS / 5 ))
[[ "$THREADS_PER_POP" -lt 1 ]] && THREADS_PER_POP=1

POPS=(AFR NFE AMR EAS SAS)
VCF_BASE="$PROJ/data/vcf"
SAMP_DIR="$PROJ/data/samples/GNOMAD_GENOME"
FREQ_DIR="$WORK/freq_genome"
COMMON_DIR="$WORK/freq_common_genome"
FINAL_PANEL="$REF/FINAL_FREQUENCIES_GNOMAD_GENOME.txt"
METADATA="$PROJ/data/panel/gnomad_hgdp_1kg_meta.tsv"
N_SAMPLES=100
MIN_AF=0.01

mkdir -p "$SAMP_DIR" "$FREQ_DIR" "$COMMON_DIR" "$REF" "$WORK/tmp"

echo "============================================================"
echo " Panel gnomAD N5 — Genoma completo"
echo " Poblaciones : ${POPS[*]}"
echo " Muestras    : $N_SAMPLES por población"
echo " Threads     : $THREADS ($THREADS_PER_POP por pop)"
echo " Panel final : $FINAL_PANEL"
echo "============================================================"

# -----------------------------------------------------------------------------
# PASO 1: Selección de muestras
# Aleatoria, high_quality=true, sin relacionados, verificadas en VCF
# Se hace una sola vez usando chr22 como referencia de IDs
# -----------------------------------------------------------------------------
echo ""
echo "[1/4] Seleccionando muestras..."

VCF_REF="$VCF_BASE/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz"
VCF_SAMPLES="$SAMP_DIR/vcf_samples.txt"

if [[ ! -s "$VCF_SAMPLES" ]]; then
  echo "  Extrayendo IDs del VCF de referencia (chr22)..."
  bcftools query -l "$VCF_REF" > "$VCF_SAMPLES"
fi
echo "  Muestras en VCF: $(wc -l < "$VCF_SAMPLES")"

for pop in "${POPS[@]}"; do
  pop_lower="${pop,,}"
  OUT_SAMP="$SAMP_DIR/${pop}_samples.txt"

  if [[ -s "$OUT_SAMP" ]]; then
    echo "  [=] $pop: ya seleccionado ($(wc -l < "$OUT_SAMP") muestras)"
    continue
  fi

  echo "  Seleccionando $pop..."

  CANDIDATES=$(awk -F'\t' -v p="$pop_lower" '
    NR>1 &&
    $160=="true" &&
    $118==p &&
    $140=="false"
    { print $1 }
  ' "$METADATA" | grep -Fxf "$VCF_SAMPLES" || true)

  N_CAND=$(echo "$CANDIDATES" | grep -c . 2>/dev/null || echo 0)
  echo "  $pop candidatos disponibles: $N_CAND"

  if [[ "$N_CAND" -lt "$N_SAMPLES" ]]; then
    echo "  [!] Solo $N_CAND para $pop — usando todos" >&2
    echo "$CANDIDATES" > "$OUT_SAMP"
  else
    echo "$CANDIDATES" | shuf | head -"$N_SAMPLES" > "$OUT_SAMP"
  fi

  echo "  [✓] $pop: $(wc -l < "$OUT_SAMP") muestras"
done

# -----------------------------------------------------------------------------
# Header del panel final
# -----------------------------------------------------------------------------
if [[ ! -s "$FINAL_PANEL" ]]; then
  {
    printf "#chrom\tposition\trsid\tA1\tA2"
    for pop in "${POPS[@]}"; do printf "\t%s" "$pop"; done
    printf "\n"
  } > "$FINAL_PANEL"
fi

# -----------------------------------------------------------------------------
# PASO 2-4: Por cromosoma — descargar, extraer, construir panel, borrar VCF
# -----------------------------------------------------------------------------
GNOMAD_BASE="https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes"

for chr in $(seq 1 22); do
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo " chr${chr}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  VCF_NAME="gnomad.genomes.v3.1.2.hgdp_tgp.chr${chr}.vcf.bgz"
  VCF_PATH="$VCF_BASE/$VCF_NAME"

  # Verificar si este cromosoma ya fue procesado (hay TSVs para todas las pops)
  all_done=true
  for pop in "${POPS[@]}"; do
    if [[ ! -s "$FREQ_DIR/${pop}_chr${chr}_af.tsv" ]]; then
      all_done=false
      break
    fi
  done

  if [[ "$all_done" == true ]]; then
    echo "  [=] chr${chr} ya procesado — saltando descarga."
  else
    # Descargar VCF si no existe (chr22 ya lo tenemos)
    if [[ ! -s "$VCF_PATH" ]]; then
      echo "  [↓] Descargando $VCF_NAME..."
      echo "  Espacio disponible: $(df -h "$VCF_BASE" | tail -1 | awk '{print $4}')"

      wget -c -q --show-progress \
        "$GNOMAD_BASE/${VCF_NAME}" \
        -O "$VCF_PATH"
      wget -c -q --show-progress \
        "$GNOMAD_BASE/${VCF_NAME}.tbi" \
        -O "${VCF_PATH}.tbi"

      echo "  [✓] Descarga completada: $(ls -lh "$VCF_PATH" | awk '{print $5}')"
    fi

    # ------------------------------------------------------------------
    # PASO 2: Extraer AF por población en paralelo (5 pops simultáneas)
    # ------------------------------------------------------------------
    echo "  [2/4] Extrayendo frecuencias alélicas..."

    extract_freq() {
      local vcf="$1"
      local pop="$2"
      local samples_file="$3"
      local out_tsv="$4"
      local threads="$5"
      local chr="$6"

      if [[ -s "$out_tsv" ]]; then
        echo "    [=] $pop chr${chr} ya procesado."
        return
      fi

      echo "    [..] $pop chr${chr}..."

      bcftools view --threads "$threads" \
          -S "$samples_file" -f PASS -m2 -M2 -v snps "$vcf" \
        | bcftools +fill-tags -- -t AF \
        | bcftools query \
            -i "AF >= $MIN_AF && AF <= $(echo "1 - $MIN_AF" | bc)" \
            -f "%CHROM\t%POS\t%REF\t%ALT\t%AF\n" \
        | awk -v maf="$MIN_AF" 'BEGIN{OFS="\t"}
            NF==5 {
              chr=$1; pos=$2; ref=$3; alt=$4; af=$5
              if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt && af+0>=maf+0) {
                if (af < 0.0001) af = 0.0001
                if (af > 0.9999) af = 0.9999
                printf "%s\t%s\t%s\t%s\t%.6f\n", chr, pos, ref, alt, af
              }
            }' > "$out_tsv"

      echo "    [✓] $pop chr${chr}: $(wc -l < "$out_tsv" | tr -d ' ') SNPs"
    }

    export -f extract_freq
    export MIN_AF

    pids=()
    for pop in "${POPS[@]}"; do
      SAMP="$SAMP_DIR/${pop}_samples.txt"
      OUT_TSV="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
      extract_freq "$VCF_PATH" "$pop" "$SAMP" "$OUT_TSV" "$THREADS_PER_POP" "$chr" &
      pids+=($!)
    done

    for pid in "${pids[@]}"; do
      wait "$pid" || { echo "[!] Error en una población chr${chr}" >&2; }
    done

    # ------------------------------------------------------------------
    # Borrar VCF de gnomAD para liberar disco
    # NO borrar chr22 porque ya lo teníamos antes
    # ------------------------------------------------------------------
    if [[ "$chr" -ne 22 ]]; then
      echo "  [🗑] Borrando $VCF_NAME para liberar disco..."
      rm -f "$VCF_PATH" "${VCF_PATH}.tbi"
      echo "  Espacio disponible ahora: $(df -h "$VCF_BASE" | tail -1 | awk '{print $4}')"
    fi
  fi

  # Verificar que todos los TSVs existen
  all_ok=true
  for pop in "${POPS[@]}"; do
    if [[ ! -s "$FREQ_DIR/${pop}_chr${chr}_af.tsv" ]]; then
      echo "  [!] Sin datos para $pop chr${chr} — saltando" >&2
      all_ok=false
    fi
  done
  [[ "$all_ok" == false ]] && continue

  # ------------------------------------------------------------------
  # PASO 3: Intersección SNPs comunes en las 5 poblaciones
  # ------------------------------------------------------------------
  echo "  [3/4] Intersección SNPs comunes..."

  id_files=()
  for pop in "${POPS[@]}"; do
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

  N_POPS=${#POPS[@]}
  COMMON_IDS="$COMMON_DIR/chr${chr}.ids"

  cat "${id_files[@]}" | sort | awk -v n="$N_POPS" '
    { if ($0==prev){c++} else { if (prev!="" && c==n) print prev; prev=$0; c=1 } }
    END{ if (prev!="" && c==n) print prev }
  ' > "$COMMON_IDS"

  N_COMMON=$(wc -l < "$COMMON_IDS" | tr -d ' ')
  echo "  [✓] SNPs comunes chr${chr}: $N_COMMON"
  [[ "$N_COMMON" -eq 0 ]] && continue

  # ------------------------------------------------------------------
  # PASO 4: Construir filas del panel
  # ------------------------------------------------------------------
  echo "  [4/4] Añadiendo chr${chr} al panel..."

  files=("$COMMON_IDS")
  for pop in "${POPS[@]}"; do files+=("$COMMON_DIR/${pop}_chr${chr}.tab"); done

  awk -v pops_str="${POPS[*]}" '
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
    }' "${files[@]}" \
    | sort -k1,1 -k2,2n \
    >> "$FINAL_PANEL"

  echo "  [✓] chr${chr} añadido."

  # Limpiar intermedios del cromosoma para no acumular
  rm -f "$COMMON_DIR/"*"_chr${chr}"*

done

# -----------------------------------------------------------------------------
# Validación final
# -----------------------------------------------------------------------------
NF_EXPECTED=$((5 + ${#POPS[@]}))
awk -v nf="$NF_EXPECTED" 'BEGIN{FS=OFS="\t"} NR==1{print;next} NF==nf{print}' \
  "$FINAL_PANEL" > "${FINAL_PANEL}.tmp" && mv "${FINAL_PANEL}.tmp" "$FINAL_PANEL"

NSNPS=$(awk 'NR>1{c++} END{print c+0}' "$FINAL_PANEL")

echo ""
echo "============================================================"
echo " ✅ Panel completado"
echo "    Archivo    : $FINAL_PANEL"
echo "    SNPs totales: $NSNPS"
echo "    Poblaciones : ${POPS[*]}"
echo "    Muestras    : $N_SAMPLES por población"
echo "    Filtro AF   : >= $MIN_AF"
echo "    Genoma      : GRCh38"
echo "============================================================"
