#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_DIR="$PROJ_DIR/data/samples"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

N_RUNS=10

# 11 poblaciones
POPULATIONS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

mkdir -p "$SAMPLES_DIR" "$OUTPUT_DIR" "$WORK_DIR/panel_N11_runs"

echo "=============================================="
echo "PANEL B: N11 (11 poblaciones)"
echo "10 corridas × muestreo SIN reemplazo"
echo "=============================================="

# Extraer listas y determinar N por población
declare -A POP_N
for pop in "${POPULATIONS[@]}"; do
    awk -v pop="$pop" '$2==pop {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_${pop}.txt"
    n=$(wc -l < "$SAMPLES_DIR/all_${pop}.txt")
    
    # Usar mínimo entre 100 y disponibles
    if [ $n -ge 100 ]; then
        POP_N[$pop]=100
    else
        POP_N[$pop]=$n
    fi
    
    echo "  $pop: $n disponibles → usaremos ${POP_N[$pop]} por corrida"
done
echo ""

# Función: seleccionar N al azar SIN reemplazo
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
    local run_dir=$3
    local output_prefix=$4
    
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    echo "  [Run$run Chr$chr]"
    
    # Procesar cada población
    for pop in "${POPULATIONS[@]}"; do
        local samples="$run_dir/samples_${pop}.txt"
        local pop_vcf="${output_prefix}_${pop}.vcf.gz"
        local pop_af="${output_prefix}_${pop}_AF.vcf.gz"
        local pop_freq="${output_prefix}_${pop}.freq"
        
        # Extraer + calcular AF
        bcftools view -S "$samples" -Oz -o "$pop_vcf" "$vcf"
        bcftools +fill-tags "$pop_vcf" -Oz -o "$pop_af" -- -t AF
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "$pop_af" > "$pop_freq"
        
        rm -f "$pop_vcf" "$pop_af"
    done
    
    # Combinar frecuencias
    paste <(cut -f1-5 "${output_prefix}_IBS.freq") \
          <(cut -f6 "${output_prefix}_IBS.freq") \
          <(cut -f6 "${output_prefix}_CEU.freq") \
          <(cut -f6 "${output_prefix}_YRI.freq") \
          <(cut -f6 "${output_prefix}_ESN.freq") \
          <(cut -f6 "${output_prefix}_GWD.freq") \
          <(cut -f6 "${output_prefix}_LWK.freq") \
          <(cut -f6 "${output_prefix}_PEL.freq") \
          <(cut -f6 "${output_prefix}_MXL.freq") \
          <(cut -f6 "${output_prefix}_PUR.freq") \
          <(cut -f6 "${output_prefix}_FIN.freq") \
          <(cut -f6 "${output_prefix}_CLM.freq") \
          > "${output_prefix}.freq"
    
    rm -f "${output_prefix}"_*.freq
}

# BUCLE PRINCIPAL
for run in {1..10}; do
    echo ""
    echo "========== CORRIDA $run/$N_RUNS =========="
    
    run_dir="$WORK_DIR/panel_N11_runs/run${run}"
    mkdir -p "$run_dir"
    
    echo "Seleccionando muestras..."
    for pop in "${POPULATIONS[@]}"; do
        select_samples "$SAMPLES_DIR/all_${pop}.txt" ${POP_N[$pop]} "$run_dir/samples_${pop}.txt"
        echo "  $pop: $(wc -l < $run_dir/samples_${pop}.txt)"
    done
    echo ""
    
    for chr in {1..22}; do
        process_chromosome $chr $run "$run_dir" "$run_dir/chr${chr}"
    done
    
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tIBS\tCEU\tYRI\tESN\tGWD\tLWK\tPEL\tMXL\tPUR\tFIN\tCLM"
        cat "$run_dir"/chr{1..22}.freq
    } > "$OUTPUT_DIR/PANEL_N11_run${run}.txt"
    
    echo "✓ Run${run}: $(wc -l < $OUTPUT_DIR/PANEL_N11_run${run}.txt) SNPs"
    rm -rf "$run_dir"
done

echo ""
echo "========== PANEL B COMPLETO =========="
ls -lh "$OUTPUT_DIR"/PANEL_N11_run*.txt
