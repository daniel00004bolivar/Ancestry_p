#!/usr/bin/env bash
# =============================================================================
# build_panel_streaming.sh
# Panel gnomAD v3.1.2 sites VCF — sin cargar todo en RAM
# Usa sort+join en disco en lugar de awk en memoria
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

POPS_ORDER=(AFR NFE AMR EAS SAS)
declare -A POP_AC=([AFR]="AC_afr" [NFE]="AC_nfe" [AMR]="AC_amr" [EAS]="AC_eas" [SAS]="AC_sas")
declare -A POP_AN=([AFR]="AN_afr" [NFE]="AN_nfe" [AMR]="AN_amr" [EAS]="AN_eas" [SAS]="AN_sas")

GNOMAD_BASE="https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes"
FINAL_PANEL="$REF/FINAL_FREQUENCIES_GNOMAD3.txt"
FREQ_DIR="$WORK/freq_raw"
COMMON_DIR="$WORK/freq_common"
VCF_DIR="$PROJ_DIR/data/vcf"
THREADS="${1:-20}"

mkdir -p "$FREQ_DIR" "$COMMON_DIR" "$VCF_DIR" "$REF"

# Header
{
  printf "#chrom\tposition\trsid\tA1\tA2"
  for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done
  printf "\n"
} > "$FINAL_PANEL"

echo "============================================================"
echo " Panel gnomAD v3.1.2 — streaming (sin carga en RAM)"
echo " Panel: $FINAL_PANEL"
echo "============================================================"

for chr in $(seq 1 22); do
  echo ""
  echo "── chr${chr} ──────────────────────────────────────"

  VCF_NAME="gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz"
  VCF_PATH="$VCF_DIR/$VCF_NAME"

  # Verificar si ya procesamos este cromosoma
  already_done=true
  for pop in "${POPS_ORDER[@]}"; do
    [[ ! -s "$FREQ_DIR/${pop}_chr${chr}_af.tsv" ]] && already_done=false && break
  done

  if [[ "$already_done" == true ]]; then
    echo "  [=] chr${chr} ya procesado — saltando descarga."
  else
    # Descargar si no existe
    if [[ ! -s "$VCF_PATH" ]]; then
      echo "  [↓] Descargando $VCF_NAME..."
      wget -c -q --show-progress \
        "$GNOMAD_BASE/$VCF_NAME" -O "$VCF_PATH"
      wget -c -q --show-progress \
        "$GNOMAD_BASE/${VCF_NAME}.tbi" -O "${VCF_PATH}.tbi"
    fi

    # Extraer AF para cada población en paralelo
    # Cada población: streaming con bcftools query, sin cargar en RAM
    echo "  [..] Extrayendo frecuencias chr${chr}..."
    pids=()
    for pop in "${POPS_ORDER[@]}"; do
      out="$FREQ_DIR/${pop}_chr${chr}_af.tsv"
      [[ -s "$out" ]] && continue
      ac="${POP_AC[$pop]}"
      an="${POP_AN[$pop]}"
      (
        bcftools view --threads 4 -f PASS -m2 -M2 -v snps "$VCF_PATH" \
          | bcftools query \
              -i "${ac} != \".\" && ${an} > 0 && ${ac}/${an} >= 0.0001" \
              -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/${ac}\t%INFO/${an}\n" \
          | awk 'BEGIN{OFS="\t"}
            NF==6 {
              chr=$1; sub(/^chr/,"",chr)
              pos=$2; ref=$3; alt=$4; ac=$5; an=$6
              if (ref~/^[ACGT]$/ && alt~/^[ACGT]$/ && ref!=alt && an>0 && ac>=0) {
                af = ac/an
                if (af<0.0001) af=0.0001
                if (af>0.9999) af=0.9999
                printf "%s\t%s\t%s\t%s\t%.6f\n", chr, pos, ref, alt, af
              }
            }' | awk '$5+0 > 0.001' | sort -k1,1 -k2,2n -k3,3 -k4,4 \
          > "$out"
        echo "  [✓] $pop chr${chr}: $(wc -l < "$out") SNPs"
      ) &
      pids+=($!)
    done
    for pid in "${pids[@]}"; do wait "$pid"; done

    # Borrar VCF para liberar disco
    if [[ "$chr" -ne 22 ]]; then
      echo "  [🗑] Borrando $VCF_NAME..."
      rm -f "$VCF_PATH" "${VCF_PATH}.tbi"
    fi
  fi

  # Verificar TSVs
  all_ok=true
  for pop in "${POPS_ORDER[@]}"; do
    [[ ! -s "$FREQ_DIR/${pop}_chr${chr}_af.tsv" ]] && { echo "  [!] Falta $pop chr${chr}"; all_ok=false; }
  done
  [[ "$all_ok" == false ]] && continue

  # ── Join por streaming usando sort+join ──────────────────
  echo "  [..] Construyendo panel chr${chr} por streaming..."

  # Crear clave chr:pos:ref:alt para cada población
  for pop in "${POPS_ORDER[@]}"; do
    awk 'BEGIN{OFS="\t"} {print $1":"$2":"$3":"$4, $5}' \
      "$FREQ_DIR/${pop}_chr${chr}_af.tsv" \
      | sort -k1,1 \
      > "$COMMON_DIR/${pop}_chr${chr}_keyed.tsv"
  done

  # Join secuencial de las 5 poblaciones
  # Empezar con AFR, ir haciendo inner join con cada población
  cp "$COMMON_DIR/AFR_chr${chr}_keyed.tsv" "$COMMON_DIR/joined_chr${chr}.tsv"

  for pop in NFE AMR EAS SAS; do
    join -t $'\t' -j 1 \
      "$COMMON_DIR/joined_chr${chr}.tsv" \
      "$COMMON_DIR/${pop}_chr${chr}_keyed.tsv" \
      > "$COMMON_DIR/joined_chr${chr}_tmp.tsv"
    mv "$COMMON_DIR/joined_chr${chr}_tmp.tsv" "$COMMON_DIR/joined_chr${chr}.tsv"
  done

  N_COMMON=$(wc -l < "$COMMON_DIR/joined_chr${chr}.tsv")
  echo "  [✓] SNPs comunes chr${chr}: $N_COMMON"
  [[ "$N_COMMON" -eq 0 ]] && continue

  # Convertir a formato panel: chr pos rsid ref alt AF_AFR AF_NFE AF_AMR AF_EAS AF_SAS
  awk 'BEGIN{OFS="\t"}
    {
      split($1, k, ":")
      chr=k[1]; pos=k[2]; ref=k[3]; alt=k[4]
      rsid="rs"chr"_"pos"_"ref"_"alt
      # campos: key AF_AFR AF_NFE AF_AMR AF_EAS AF_SAS
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        chr, pos, rsid, ref, alt, $2, $3, $4, $5, $6
    }' "$COMMON_DIR/joined_chr${chr}.tsv" \
    | sort -k1,1n -k2,2n \
    >> "$FINAL_PANEL"

  echo "  [✓] chr${chr} añadido al panel."

  # Limpiar intermedios del cromosoma
  rm -f "$COMMON_DIR/"*"_chr${chr}"*

done

# Validación final
NF_EXPECTED=10
awk -v nf="$NF_EXPECTED" 'BEGIN{FS=OFS="\t"} NR==1{print;next} NF==nf{print}' \
  "$FINAL_PANEL" > "${FINAL_PANEL}.tmp" && mv "${FINAL_PANEL}.tmp" "$FINAL_PANEL"

NSNPS=$(awk 'NR>1{c++} END{print c+0}' "$FINAL_PANEL")

echo ""
echo "============================================================"
echo " ✅ Panel completado"
echo "    Archivo : $FINAL_PANEL"
echo "    SNPs    : $NSNPS"
echo "    Genoma  : GRCh38"
echo "============================================================"
