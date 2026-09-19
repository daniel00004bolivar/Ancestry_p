#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   bash bin/build_final_from_af_tsv.sh <POPS_FILE> <CHR_LIST>
#
# Lee work/freq_raw/<POP>_chr<CHR>_af.tsv (CHROM POS REF ALT AF)
# Construye refpanels/FINAL_FREQUENCIES.txt con:
# #chrom position rsid A1 A2 POP1 POP2...

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
read -r -a CHRS <<< "$CHR_LIST"

FINAL="$REF/FINAL_FREQUENCIES.txt"
mkdir -p "$WORK/freq_common" "$REF"

# Header
{
  printf "#chrom\tposition\trsid\tA1\tA2"
  for p in "${POPS[@]}"; do printf "\t%s" "$p"; done
  printf "\n"
} > "$FINAL"

for chr in "${CHRS[@]}"; do
  echo "Procesando chr${chr}…"

  # 1) Construye IDS por pop, luego intersección
  id_files=()
  for pop in "${POPS[@]}"; do
    in="$WORK/freq_raw/${pop}_chr${chr}_af.tsv"
    [[ -s "$in" ]] || { echo "[X] Falta $in" >&2; exit 1; }

    tab="$WORK/freq_common/${pop}_chr${chr}.tab"
    ids="$WORK/freq_common/${pop}_chr${chr}.ids"

    # key, chr, pos, ref, alt, af
    awk 'BEGIN{OFS="\t"}
      NF>=5 {
        chr=$1; pos=$2; ref=$3; alt=$4; af=$5
        if (ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/ && ref!=alt) {
          key=chr ":" pos ":" ref ":" alt
          print key, chr, pos, ref, alt, af
        }
      }' "$in" > "$tab"

    cut -f1 "$tab" | sort -u > "$ids"
    id_files+=("$ids")
  done

  N=${#POPS[@]}
  common_ids="$WORK/freq_common/chr${chr}.ids"
  cat "${id_files[@]}" | sort | awk -v n="$N" '
    { if ($0==p){c++} else { if (p!="" && c==n) print p; p=$0; c=1 } }
    END{ if (p!="" && c==n) print p }' > "$common_ids"

  if [[ ! -s "$common_ids" ]]; then
    echo "  [!] Sin SNPs comunes en chr${chr}"
    continue
  fi

  # 2) Construye filas para ese chr
  files=("$common_ids")
  for pop in "${POPS[@]}"; do files+=("$WORK/freq_common/${pop}_chr${chr}.tab"); done

  awk -v pops="$(IFS=' '; echo "${POPS[*]}")" '
    BEGIN{OFS="\t"; split(pops,P," "); np=length(P)}
    ARGIND==1{ids[$1]=1; next}
    ARGIND>1{
      key=$1; chr=$2; pos=$3; ref=$4; alt=$5; af=$6
      idx=ARGIND-1; pop=P[idx]
      freq[pop,key]=af
      if(!(key in meta)) meta[key]=chr OFS pos OFS ("rs" chr "_" pos "_" ref "_" alt) OFS ref OFS alt
      next
    }
    END{
      for(k in ids){
        if(!(k in meta)) continue
        printf "%s", meta[k]
        for(i=1;i<=np;i++){
          pop=P[i]
          val = ((pop,k) in freq ? freq[pop,k] : "")
          # clamp 0/1 para evitar problemas numéricos:
          if (val=="" ) val=0.0001
          if (val<0.0001) val=0.0001
          if (val>0.9999) val=0.9999
          printf OFS "%.6f", val
        }
        printf "\n"
      }
    }' "${files[@]}" \
    | sort -k1,1n -k2,2n \
    >> "$FINAL"

  echo "  [✓] chr${chr} → SNPs comunes: $(wc -l < "$common_ids" | tr -d ' ')"
done

# Limpia líneas vacías / asegura columnas correctas
NF_expected=$((5 + ${#POPS[@]}))
awk -v nf="$NF_expected" 'BEGIN{FS=OFS="\t"} NR==1{print; next} NF==nf {print}' "$FINAL" > "$FINAL.tmp" && mv "$FINAL.tmp" "$FINAL"

NSNPS=$(awk 'NR>1{c++} END{print c+0}' "$FINAL")
echo "✅ Panel final creado: $FINAL"
echo "   - SNPs: $NSNPS"

