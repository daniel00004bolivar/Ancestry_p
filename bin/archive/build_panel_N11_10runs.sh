#!/bin/bash
# ==============================================================================
# PANEL B: N11 (11 poblaciones específicas)
# 10 corridas × 100 individuos por población
# Poblaciones: IBS, CEU, YRI, ESN, GWD, LWK, PEL, MXL, PUR, FIN, CLM
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

# 11 poblaciones
POPULATIONS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

mkdir -p "$SAMPLES_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR/panel_N11_runs"

echo "=============================================="
echo "PANEL B: N11 (11 poblaciones)"
echo "10 corridas × 100 muestras por población"
echo "=============================================="
echo ""

# Extraer listas de muestras por población
echo "Extrayendo listas de muestras..."
for pop in "${POPULATIONS[@]}"; do
    awk -v pop="$pop" '$2==pop {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_${pop}.txt"
    n=$(wc -l < "$SAMPLES_DIR/all_${pop}.txt")
    echo "  $pop: $n muestras disponibles"
    
    if [ "$n" -lt 100 ]; then
        echo "    ⚠️  ADVERTENCIA: Menos de 100 muestras. Habrá repeticiones."
    fi
done
echo ""

# Función para seleccionar muestras aleatorias CON reemplazo
select_samples_with_replacement() {
    local input_file=$1
    local n_samples=$2
    local seed=$3
    local output_file=$4
    
    shuf -r -n "$n_samples" --random-source=<(yes $seed) "$input_file" > "$output_file"
}

# Función para procesar un cromosoma
process_chromosome() {
    local chr=$1
    local run=$2
    local run_dir=$3
    local output_prefix=$4
    
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    if [ ! -f "$vcf" ]; then
        echo "ERROR: No se encuentra $vcf"
        return 1
    fi
    
    echo "  [Run$run Chr$chr] Procesando 11 poblaciones..."
    
    # Procesar cada población
    declare -a freq_files
    for pop in "${POPULATIONS[@]}"; do
        local samples="$run_dir/samples_${pop}.txt"
        local pop_vcf="${output_prefix}_${pop}.vcf.gz"
        local pop_af="${output_prefix}_${pop}_AF.vcf.gz"
        local pop_freq="${output_prefix}_${pop}.freq"
        
        # Extraer muestras y calcular AF
        bcftools view -S "$samples" -Oz -o "$pop_vcf" "$vcf"
        bcftools +fill-tags "$pop_vcf" -Oz -o "$pop_af" -- -t AF
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "$pop_af" > "$pop_freq"
        
        freq_files+=("$pop_freq")
        
        # Limpiar VCFs temporales
        rm -f "$pop_vcf" "$pop_af"
    done
    
    # Combinar todas las frecuencias
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
    
    # Limpiar archivos de frecuencias temporales
    rm -f "${output_prefix}"_*.freq
    
    echo "  [Run$run Chr$chr] ✓ Completado"
}

# BUCLE PRINCIPAL: 10 corridas
for run in {1..10}; do
    echo ""
    echo "=============================================="
    echo "CORRIDA $run de $N_RUNS"
    echo "=============================================="
    
    run_dir="$WORK_DIR/panel_N11_runs/run${run}"
    mkdir -p "$run_dir"
    
    # Seleccionar muestras aleatorias para cada población
    echo "Seleccionando muestras aleatorias..."
    for pop in "${POPULATIONS[@]}"; do
        select_samples_with_replacement \
            "$SAMPLES_DIR/all_${pop}.txt" \
            $N_SAMPLES \
            $run \
            "$run_dir/samples_${pop}.txt"
        
        n=$(wc -l < "$run_dir/samples_${pop}.txt")
        echo "  $pop: $n muestras"
    done
    echo ""
    
    # Procesar chr1-22
    echo "Procesando cromosomas 1-22..."
    for chr in {1..22}; do
        process_chromosome $chr $run "$run_dir" "$run_dir/chr${chr}"
    done
    
    # Concatenar todos los cromosomas
    echo "Concatenando cromosomas..."
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tIBS\tCEU\tYRI\tESN\tGWD\tLWK\tPEL\tMXL\tPUR\tFIN\tCLM"
        cat "$run_dir"/chr{1..22}.freq
    } > "$OUTPUT_DIR/PANEL_N11_run${run}.txt"
    
    # Verificar
    n_snps=$(wc -l < "$OUTPUT_DIR/PANEL_N11_run${run}.txt")
    echo "✓ Panel run${run} completo: $n_snps SNPs"
    
    # Limpiar archivos intermedios
    rm -rf "$run_dir"
done

echo ""
echo "=============================================="
echo "PANEL B COMPLETADO"
echo "=============================================="
echo ""
echo "10 paneles creados:"
ls -lh "$OUTPUT_DIR"/PANEL_N11_run*.txt
echo ""
echo "Siguiente paso: Ejecutar iAdmix con estos paneles"
