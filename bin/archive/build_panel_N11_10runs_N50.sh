#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_DIR="$PROJ_DIR/data/samples"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

N_RUNS=10
N_SAMPLES=50

POPULATIONS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

mkdir -p "$SAMPLES_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR/panel_N11_runs"

echo "=============================================="
echo "PANEL N11: 11 poblaciones"
echo "10 corridas × 50 muestras por población"
echo "=============================================="
echo ""

# Extraer listas de muestras
echo "Extrayendo listas de muestras..."
for pop in "${POPULATIONS[@]}"; do
    awk -v pop="$pop" '$2==pop {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_${pop}.txt"
    n=$(wc -l < "$SAMPLES_DIR/all_${pop}.txt")
    echo "  $pop: $n muestras disponibles"
done
echo ""

# Función para seleccionar muestras SIN reemplazo
select_samples() {
    local input_file=$1
    local n_samples=$2
    local seed=$3
    local output_file=$4
    
    shuf -n "$n_samples" --random-source=<(yes $seed) "$input_file" > "$output_file"
}

# Función para procesar un cromosoma
process_chromosome() {
    local chr=$1
    local run=$2
    local run_dir=$3
    local output_prefix=$4
    
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    echo "  [Run$run Chr$chr]"
    
    for pop in "${POPULATIONS[@]}"; do
        samples="$run_dir/samples_${pop}.txt"
        bcftools view -S "$samples" -Oz -o "${output_prefix}_${pop}.vcf.gz" "$vcf"
        bcftools +fill-tags "${output_prefix}_${pop}.vcf.gz" -Oz -o "${output_prefix}_${pop}_AF.vcf.gz" -- -t AF
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' "${output_prefix}_${pop}_AF.vcf.gz" > "${output_prefix}_${pop}.freq"
    done
    
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
    
    rm -f "${output_prefix}"_{IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM}{,_AF}.vcf.gz "${output_prefix}"_{IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM}.freq
}

# BUCLE PRINCIPAL
for run in {1..10}; do
    echo ""
    echo "========== CORRIDA $run/$N_RUNS =========="
    
    run_dir="$WORK_DIR/panel_N11_runs/run${run}"
    mkdir -p "$run_dir"
    
    echo "Seleccionando muestras (n=50, SIN reemplazo)..."
    for pop in "${POPULATIONS[@]}"; do
        select_samples \
            "$SAMPLES_DIR/all_${pop}.txt" \
            $N_SAMPLES \
            $run \
            "$run_dir/samples_${pop}.txt"
        
        n=$(wc -l < "$run_dir/samples_${pop}.txt")
        echo "  $pop: $n muestras"
    done
    echo ""
    
    echo "Procesando cromosomas..."
    for chr in {1..22}; do
        process_chromosome $chr $run "$run_dir" "$run_dir/chr${chr}"
    done
    
    echo "Concatenando cromosomas..."
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tIBS\tCEU\tYRI\tESN\tGWD\tLWK\tPEL\tMXL\tPUR\tFIN\tCLM"
        for chr in {1..22}; do
            tail -n +2 "$run_dir/chr${chr}.freq"
        done
    } > "$OUTPUT_DIR/PANEL_N11_run${run}.txt"
    
    rm -rf "$run_dir/chr"*.freq "$run_dir/chr"*.vcf.gz
    
    echo "✓ Panel N11 run$run completado"
done

echo ""
echo "=============================================="
echo "✓ 10 PANELES N11 COMPLETADOS"
echo "=============================================="
ls -lh "$OUTPUT_DIR"/PANEL_N11_run*.txt
