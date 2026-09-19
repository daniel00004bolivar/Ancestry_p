#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   bash bin/01_make_af_if_needed.sh <POPS_FILE> <CHR_LIST>
# Ej: bash bin/01_make_af_if_needed.sh config/pops_n5.txt "21 22"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

if [[ $# -ne 2 ]]; then
  echo "[X] Uso: $0 <POPS_FILE> <CHR_LIST>" >&2
  exit 1
fi

POPS_FILE="$1"
CHR_LIST="$2"

[[ -s "$POPS_FILE" ]] || { echo "[X] No existe/está vacío: $POPS_FILE" >&2; exit 1; }

mapfile -t POPS < <(awk 'NF{gsub("\r",""); gsub(/^[ \t]+|[ \t]+$/,""); print toupper($0)}' "$POPS_FILE")

# Normaliza CHR_LIST a array
read -r -a CHRS <<< "$CHR_LIST"

echo "Generando AF TSV para POPS: ${POPS[*]}"
echo "Cromosomas: ${CHRS[*]}"

for chr in "${CHRS[@]}"; do
  # Busca VCF del cromosoma
  vcf=""
  for pat in \
    "ALL.chr${chr}.vcf.gz" \
    "ALL.chr${chr}.*.vcf.gz" \
    "chr${chr}*.vcf.gz" \
    "*chr${chr}*.vcf.gz"
  do
    f=$(ls -1 "$VCF_DIR"/$pat 2>/dev/null | head -n1 || true)
    if [[ -n "${f:-}" ]]; then vcf="$f"; break; fi
  done

  if [[ -z "$vcf" ]]; then
    echo "[X] No VCF chr${chr} en $VCF_DIR" >&2
    exit 1
  fi

  # Verifica índice .tbi
  if [[ ! -s "${vcf}.tbi" ]]; then
    echo "[X] Falta índice .tbi para: $vcf" >&2
    echo "    Ejecuta: tabix -p vcf \"$vcf\"" >&2
    exit 1
  fi

  for pop in "${POPS[@]}"; do
    samp="$SAMP_DIR/${pop}_samples.txt"
    out="$WORK/freq_raw/${pop}_chr${chr}_af.tsv"

    if [[ ! -s "$samp" ]]; then
      echo "[!] No hay muestras para $pop, se omite."
      continue
    fi

    if [[ -s "$out" ]]; then
      echo "[=] $pop chr${chr} ya procesado."
      continue
    fi

    echo "[..] $pop chr${chr} → AF TSV"
    # Calcula AF con bcftools:
    #  - filtra a bialélicos SNPs
    #  - rellena TAG AF
    #  - saca: CHROM POS REF ALT AF
    "$BCFTOOLS" view -S "$samp" -m2 -M2 -v snps -Ou "$vcf" \
      | "$BCFTOOLS" +fill-tags -Ou -- -t AF \
      | "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' \
      > "$out"

    if [[ ! -s "$out" ]]; then
      echo "[X] Error: no se generó $out" >&2
      exit 1
    fi
  done
done

echo "✅ AF TSV listos en: $WORK/freq_raw/"

