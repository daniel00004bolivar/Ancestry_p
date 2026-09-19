#!/bin/bash
# ==============================================================================
# PANEL A: 3 SUPERPOBLACIONES (AFR, EUR, AMR)
# 10 corridas × 100 individuos por superpoblación
# Muestreo CON reemplazo, selección aleatoria
# ==============================================================================

set -euo pipefail

# Configuración
PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_DIR="$PROJ_DIR/data/samples"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

N_RUNS=10
N_SAMPLES=100

mkdir -p "$SAMPLES_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR/panel_3way_runs"

echo "=============================================="
echo "PANEL A: 3 SUPERPOBLACIONES"
echo "10 corridas × 100 muestras (AFR, EUR, AMR)"
echo "=============================================="
echo ""

# Extraer listas de muestras por superpoblación
echo "Extrayendo listas de muestras..."
awk '$3=="AFR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_AFR.txt"
awk '$3=="EUR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_EUR.txt"
awk '$3=="AMR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_AMR.txt"

n_afr=$(wc -l < "$SAMPLES_DIR/all_AFR.txt")
n_eur=$(wc -l < "$SAMPLES_DIR/all_EUR.txt")
n_amr=$(wc -l < "$SAMPLES_DIR/all_AMR.txt")

echo "  AFR: $n_afr muestras disponibles"
echo "  EUR: $n_eur muestras disponibles"
echo "  AMR: $n_amr muestras disponibles"
echo ""

# Función para seleccionar muestras aleatorias CON reemplazo
select_samples_with_replacement() {
    local input_file=$1
    local n_samples=$2
    local seed=$3
    local output_file=$4
    
    # Usar shuf con reemplazo (-r) y semilla para reproducibilidad
    shuf -r -n "$n_samples" --random-source=<(yes $seed) "$input_file" > "$output_file"
}

# Función para procesar un cromosoma
process_chromosome() {
    local chr=$1
    local run=$2
    local samples_afr=$3
    local samples_eur=$4
    local samples_amr=$5
    local output_prefix=$6
    
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    if [ ! -f "$vcf" ]; then
        echo "ERROR: No se encuentra $vcf"
        return 1
    fi
    
    echo "  [Run$run Chr$chr] Procesando..."
    
    # Extraer muestras AFR
    bcftools view -S "$samples_afr" -Oz -o "${output_prefix}_AFR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_AFR.vcf.gz" -Oz -o "${output_prefix}_AFR_AF.vcf.gz" -- -t AF
    
    # Extraer muestras EUR
    bcftools view -S "$samples_eur" -Oz -o "${output_prefix}_EUR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_EUR.vcf.gz" -Oz -o "${output_prefix}_EUR_AF.vcf.gz" -- -t AF
    
    # Extraer muestras AMR
    bcftools view -S "$samples_amr" -Oz -o "${output_prefix}_AMR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_AMR.vcf.gz" -Oz -o "${output_prefix}_AMR_AF.vcf.gz" -- -t AF
    
    # Combinar frecuencias en formato iAdmix
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_AFR_AF.vcf.gz" > "${output_prefix}_AFR.freq"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_EUR_AF.vcf.gz" > "${output_prefix}_EUR.freq"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_AMR_AF.vcf.gz" > "${output_prefix}_AMR.freq"
    
    # Merge en una sola tabla
    paste <(cut -f1-5 "${output_prefix}_AFR.freq") \
          <(cut -f6 "${output_prefix}_AFR.freq") \
          <(cut -f6 "${output_prefix}_EUR.freq") \
          <(cut -f6 "${output_prefix}_AMR.freq") \
          > "${output_prefix}.freq"
    
    # Limpiar archivos temporales
    rm -f "${output_prefix}"_{AFR,EUR,AMR}{,_AF}.vcf.gz "${output_prefix}"_{AFR,EUR,AMR}.freq
    
    echo "  [Run$run Chr$chr] ✓ Completado"
}

# BUCLE PRINCIPAL: 10 corridas
for run in {1..10}; do
    echo ""
    echo "=============================================="
    echo "CORRIDA $run de $N_RUNS"
    echo "=============================================="
    
    run_dir="$WORK_DIR/panel_3way_runs/run${run}"
    mkdir -p "$run_dir"
    
    # Seleccionar muestras aleatorias (semilla = número de corrida)
    echo "Seleccionando muestras aleatorias..."
    select_samples_with_replacement "$SAMPLES_DIR/all_AFR.txt" $N_SAMPLES $run "$run_dir/samples_AFR.txt"
    select_samples_with_replacement "$SAMPLES_DIR/all_EUR.txt" $N_SAMPLES $run "$run_dir/samples_EUR.txt"
    select_samples_with_replacement "$SAMPLES_DIR/all_AMR.txt" $N_SAMPLES $run "$run_dir/samples_AMR.txt"
    
    echo "  AFR: $(wc -l < "$run_dir/samples_AFR.txt") muestras"
    echo "  EUR: $(wc -l < "$run_dir/samples_EUR.txt") muestras"
    echo "  AMR: $(wc -l < "$run_dir/samples_AMR.txt") muestras"
    echo ""
    
    # Procesar chr1-22
    echo "Procesando cromosomas 1-22..."
    for chr in {1..22}; do
        process_chromosome $chr $run \
            "$run_dir/samples_AFR.txt" \
            "$run_dir/samples_EUR.txt" \
            "$run_dir/samples_AMR.txt" \
            "$run_dir/chr${chr}"
    done
    
    # Concatenar todos los cromosomas
    echo "Concatenando cromosomas..."
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tAFR\tEUR\tAMR"
        cat "$run_dir"/chr{1..22}.freq
    } > "$OUTPUT_DIR/PANEL_3WAY_run${run}.txt"
    
    # Verificar
    n_snps=$(wc -l < "$OUTPUT_DIR/PANEL_3WAY_run${run}.txt")
    echo "✓ Panel run${run} completo: $n_snps SNPs"
    
    # Limpiar archivos intermedios
    rm -rf "$run_dir"
done

echo ""
echo "=============================================="
echo "PANEL A COMPLETADO"
echo "=============================================="
echo ""
echo "10 paneles creados:"
ls -lh "$OUTPUT_DIR"/PANEL_3WAY_run*.txt
echo ""
echo "Siguiente paso: Ejecutar iAdmix con estos paneles"
