#!/bin/bash
# 06_build_5pop_panel.sh — Construir panel de 5 superpoblaciones desde 1KGP
#
# Usa tabix para consultar el VCF de 1000 Genomes Phase 3 de forma remota,
# extrayendo solo las posiciones del panel HapMap3 (249k SNPs).
# Agrega AMR y SAS a las 3 superpoblaciones ya disponibles (AFR, EAS, EUR).
#
# Requisitos: tabix (htslib), internet, ~10 minutos
#
# METHODOLOGICAL WARNING: this panel feeds NGSadmix/PCAngsd (bin/08*,
# bin/09*), which represent each pool (POOL1/POOL2/HOSPITAL, ~50 real
# individuals each) as ONE diploid pseudo-individual -- a structural
# limitation of the Beagle format, not configurable. See
# bin/README_VIGENTE.md and bin/16_pool_replicates.sh.
#
# Usage: bash bin/06_build_5pop_panel.sh
set -euo pipefail
echo "[NOTE] Panel for NGSadmix/PCAngsd: pools = diploid pseudo-individual (see bin/README_VIGENTE.md)" >&2

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_PATTERN="${BASE_URL}/ALL.chr{CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
HAPMAP_POS="/tmp/hapmap3_pos.txt"
OUT_DIR="data/iadmix/DATA"
TMP_DIR="/tmp/1kgp_extract"
FINAL_PANEL="${OUT_DIR}/hapmap3_5pops.freqs"

mkdir -p "$TMP_DIR"

# Verificar que exista el archivo de posiciones HapMap3
if [ ! -f "$HAPMAP_POS" ]; then
    echo "Creando archivo de posiciones HapMap3..."
    awk 'NR>1{chr=$1; sub(/^chr/,"",chr); print chr"\t"$2}' \
        "${OUT_DIR}/hapmap3.8populations.hg19.freqs" > "$HAPMAP_POS"
fi

echo "Posiciones HapMap3: $(wc -l < "$HAPMAP_POS") SNPs"
echo ""

# Extraer datos de 1KGP por cromosoma (paralelo, máx 4 simultáneos)
extract_chr() {
    local CHR=$1
    local URL="${VCF_PATTERN/\{CHR\}/$CHR}"
    local REGIONS_FILE="$TMP_DIR/regions_chr${CHR}.txt"
    local OUT_FILE="$TMP_DIR/chr${CHR}_1kgp.txt"

    # Crear archivo de regiones para este cromosoma (formato: chr start end)
    awk -v c="$CHR" '$1==c {print c"\t"$2-1"\t"$2}' "$HAPMAP_POS" > "$REGIONS_FILE"

    n_pos=$(wc -l < "$REGIONS_FILE")
    if [ "$n_pos" -eq 0 ]; then
        echo "  chr${CHR}: sin posiciones — omitiendo"
        touch "$OUT_FILE"
        return
    fi

    echo "  chr${CHR}: consultando $n_pos posiciones..."
    tabix -R "$REGIONS_FILE" "$URL" 2>/dev/null | \
    awk 'BEGIN{OFS="\t"}
         /^#/ {next}
         {
            # Campos VCF: CHROM POS ID REF ALT QUAL FILTER INFO
            chr=$1; pos=$2; ref=$4; alt=$5; info=$8

            # Solo SNPs bialélicos (REF y ALT de 1 base)
            if (length(ref)!=1 || length(alt)!=1) next

            # Extraer frecuencias del campo INFO
            amr_af="."; afr_af="."; eas_af="."; eur_af="."; sas_af="."
            n=split(info, fields, ";")
            for(i=1; i<=n; i++) {
                if (fields[i] ~ /^AMR_AF=/) { amr_af=substr(fields[i],8) }
                if (fields[i] ~ /^AFR_AF=/) { afr_af=substr(fields[i],8) }
                if (fields[i] ~ /^EAS_AF=/) { eas_af=substr(fields[i],8) }
                if (fields[i] ~ /^EUR_AF=/) { eur_af=substr(fields[i],8) }
                if (fields[i] ~ /^SAS_AF=/) { sas_af=substr(fields[i],8) }
            }

            # La frecuencia en 1KGP es del alelo ALT
            # Estandarizar: freq del alelo lex-menor
            if (ref <= alt) {
                allele1=ref; allele2=alt
                # ref es lex-menor: freq_a1 = 1 - AFR_AF (ya que AFR_AF es de ALT)
                f_amr=(amr_af=="."?".":1-amr_af+0)
                f_afr=(afr_af=="."?".":1-afr_af+0)
                f_eas=(eas_af=="."?".":1-eas_af+0)
                f_eur=(eur_af=="."?".":1-eur_af+0)
                f_sas=(sas_af=="."?".":1-sas_af+0)
            } else {
                allele1=alt; allele2=ref
                # alt es lex-menor: freq_a1 = AFR_AF (frecuencia del alelo ALT)
                f_amr=amr_af; f_afr=afr_af; f_eas=eas_af; f_eur=eur_af; f_sas=sas_af
            }

            print chr, pos, allele1, allele2, f_afr, f_eas, f_eur, f_amr, f_sas
         }' > "$OUT_FILE"

    n_out=$(wc -l < "$OUT_FILE")
    echo "  chr${CHR}: $n_out SNPs extraídos"
}

export -f extract_chr
export HAPMAP_POS TMP_DIR VCF_PATTERN

echo "Extrayendo frecuencias 1KGP (tabix remoto, paralelo)..."
echo ""

# Procesar cromosomas en grupos de 4
for batch_start in 1 5 9 13 17 21; do
    batch_end=$((batch_start + 3))
    [ $batch_end -gt 22 ] && batch_end=22
    for CHR in $(seq $batch_start $batch_end); do
        extract_chr "$CHR" &
    done
    wait
done

echo ""
echo "Combinando cromosomas..."

# Combinar todos los cromosomas
echo "#chr	pos	allele1	allele2	AFR_1KGP	EAS_1KGP	EUR_1KGP	AMR_1KGP	SAS_1KGP" > "$FINAL_PANEL"
for CHR in $(seq 1 22); do
    f="$TMP_DIR/chr${CHR}_1kgp.txt"
    [ -f "$f" ] && cat "$f"
done | sort -k1,1n -k2,2n >> "$FINAL_PANEL"

total=$(wc -l < "$FINAL_PANEL")
echo "Panel combinado: $((total - 1)) SNPs → ${FINAL_PANEL}"
echo ""
echo "Listo. Ejecutar ahora:"
echo "  Rscript bin/05c_continental_ancestry_v3.R"
