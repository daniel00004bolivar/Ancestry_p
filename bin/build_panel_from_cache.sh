#!/bin/bash
set -euo pipefail
PROJ_DIR="/home/laboratorio/Documentos/Ancestria"
FINAL_PANEL="$PROJ_DIR/refpanels/FINAL_FREQUENCIES_GNOMAD3.txt"
COMMON_DIR="$PROJ_DIR/work/gnomad3_build/freq_common"
LIFTOVER_DIR="$PROJ_DIR/work/gnomad3_build/liftover"
LIFTOVER_DONE="$LIFTOVER_DIR/pos_map_37to38.tsv"
POPS_ORDER=(AFR EUR AMR EAS SAS)
CHRS=(1 2 3 4 5 6 7)
{ printf "#chrom\tposition\trsid\tA1\tA2"; for pop in "${POPS_ORDER[@]}"; do printf "\t%s" "$pop"; done; printf "\n"; } > "$FINAL_PANEL"
for chr in "${CHRS[@]}"; do
  echo "chr${chr}..."
  COMMON_IDS="$COMMON_DIR/chr${chr}.ids"
  POS_MAP="$LIFTOVER_DIR/chr${chr}_38to37.tsv"
  [[ ! -s "$COMMON_IDS" ]] && continue
  files=("$COMMON_IDS")
  for pop in "${POPS_ORDER[@]}"; do files+=("$COMMON_DIR/${pop}_chr${chr}.tab"); done
  awk -v pops_str="$(IFS=' '; echo "${POPS_ORDER[*]}")" -v pos_map_file="$POS_MAP" 'BEGIN{OFS="\t";split(pops_str,P," ");np=length(P);while((getline line<pos_map_file)>0){split(line,a,"\t");map38to37[a[1]]=a[2]}}ARGIND==1{ids[$1]=1;next}ARGIND>1{key=$1;chr=$2;pos=$3;ref=$4;alt=$5;af=$6;idx=ARGIND-1;pop=P[idx];freq[pop,key]=af;if(!(key in meta)){k38=chr":"pos;pos37=(k38 in map38to37)?map38to37[k38]:k38;rsid="rs"chr"_"pos"_"ref"_"alt;meta[key]=chr OFS pos OFS rsid OFS ref OFS alt}next}END{for(k in ids){if(!(k in meta))continue;printf "%s",meta[k];for(i=1;i<=np;i++){pop=P[i];val=((pop SUBSEP k)in freq)?freq[pop,k]:"0.000100";printf OFS"%.6f",val+0}printf "\n"}}' "${files[@]}" | sort -k1,1n -k2,2n >> "$FINAL_PANEL"
  echo "[✓] chr${chr}"
done
echo "LISTO: $FINAL_PANEL"
wc -l "$FINAL_PANEL"
