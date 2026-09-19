#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_DIR="$PROJ_DIR/data/samples"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

N_RUNS=10
N_SAMPLES=100

mkdir -p "$SAMPLES_DIR" "$OUTPUT_DIR" "$WORK_DIR/panel_3way_runs"

echo "=============================================="
echo "PANEL A: 3 SUPERPOBLACIONES"
echo "10 corridas × 100 muestras (AFR, EUR, AMR)"
echo "=============================================="

# Extraer listas
awk '$3=="AFR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_AFR.txt"
awk '$3=="EUR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_EUR.txt"
awk '$3=="AMR" {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_AMR.txt"

echo "Muestras disponibles:"
echo "  AFR: $(wc -l < $SAMPLES_DIR/all_AFR.txt)"
echo "  EUR: $(wc -l < $SAMPLES_DIR/all_EUR.txt)"
echo "  AMR: $(wc -l < $SAMPLES_DIR/all_AMR.txt)"
echo ""

# Función simple: seleccionar N al azar (sin reemplazo)
select_samples() {
    local input_file=$1
    local n_samples=$2
    local output_file=$3
    
    shuf -n "$n_samples" "$input_file" > "$output_file"
}

# Función para procesar cromosoma
process_chromosome() {
    local chr=$1
    local run=$2
    local samples_afr=$3
    local samples_eur=$4
    local samples_amr=$5
    local output_prefix=$6
    
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    echo "  [Run$run Chr$chr]"
    
    # AFR
    bcftools view -S "$samples_afr" -Oz -o "${output_prefix}_AFR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_AFR.vcf.gz" -Oz -o "${output_prefix}_AFR_AF.vcf.gz" -- -t AF
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_AFR_AF.vcf.gz" > "${output_prefix}_AFR.freq"
    
    # EUR
    bcftools view -S "$samples_eur" -Oz -o "${output_prefix}_EUR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_EUR.vcf.gz" -Oz -o "${output_prefix}_EUR_AF.vcf.gz" -- -t AF
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_EUR_AF.vcf.gz" > "${output_prefix}_EUR.freq"
    
    # AMR
    bcftools view -S "$samples_amr" -Oz -o "${output_prefix}_AMR.vcf.gz" "$vcf"
    bcftools +fill-tags "${output_prefix}_AMR.vcf.gz" -Oz -o "${output_prefix}_AMR_AF.vcf.gz" -- -t AF
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_AMR_AF.vcf.gz" > "${output_prefix}_AMR.freq"
    
    # Combinar
    paste <(cut -f1-5 "${output_prefix}_AFR.freq") \
          <(cut -f6 "${output_prefix}_AFR.freq") \
          <(cut -f6 "${output_prefix}_EUR.freq") \
          <(cut -f6 "${output_prefix}_AMR.freq") \
          > "${output_prefix}.freq"
    
    rm -f "${output_prefix}"_{AFR,EUR,AMR}{,_AF}.vcf.gz "${output_prefix}"_{AFR,EUR,AMR}.freq
}

# BUCLE PRINCIPAL
for run in {1..10}; do
    echo ""
    echo "========== CORRIDA $run/$N_RUNS =========="
    
    run_dir="$WORK_DIR/panel_3way_runs/run${run}"
    mkdir -p "$run_dir"
    
    # Seleccionar 100 muestras al azar (diferentes en cada corrida)
    select_samples "$SAMPLES_DIR/all_AFR.txt" $N_SAMPLES "$run_dir/samples_AFR.txt"
    select_samples "$SAMPLES_DIR/all_EUR.txt" $N_SAMPLES "$run_dir/samples_EUR.txt"
    select_samples "$SAMPLES_DIR/all_AMR.txt" $N_SAMPLES "$run_dir/samples_AMR.txt"
    
    echo "Muestras seleccionadas: AFR=$(wc -l < $run_dir/samples_AFR.txt) EUR=$(wc -l < $run_dir/samples_EUR.txt) AMR=$(wc -l < $run_dir/samples_AMR.txt)"
    
    for chr in {1..22}; do
        process_chromosome $chr $run \
            "$run_dir/samples_AFR.txt" \
            "$run_dir/samples_EUR.txt" \
            "$run_dir/samples_AMR.txt" \
            "$run_dir/chr${chr}"
    done
    
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tAFR\tEUR\tAMR"
        cat "$run_dir"/chr{1..22}.freq
    } > "$OUTPUT_DIR/PANEL_3WAY_run${run}.txt"
    
    echo "✓ Run${run}: $(wc -l < $OUTPUT_DIR/PANEL_3WAY_run${run}.txt) SNPs"
    rm -rf "$run_dir"
done

echo ""
echo "========== PANEL A COMPLETO =========="
ls -lh "$OUTPUT_DIR"/PANEL_3WAY_run*.txt
