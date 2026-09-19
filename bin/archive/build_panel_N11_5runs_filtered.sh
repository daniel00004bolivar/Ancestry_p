#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_DIR="$PROJ_DIR/data/samples"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

N_RUNS=5
N_SAMPLES=50
PARALLEL_CHR=6

POPULATIONS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

mkdir -p "$SAMPLES_DIR" "$OUTPUT_DIR" "$WORK_DIR/panel_N11_runs"

echo "=============================================="
echo "PANEL N11 FILTRADO: 11 poblaciones"
echo "5 corridas × 50 muestras × MAF > 0.1"
echo "Paralelizado: 6 cromosomas simultáneos"
echo "=============================================="
echo ""

# Extraer listas
echo "Extrayendo listas de muestras..."
for pop in "${POPULATIONS[@]}"; do
    awk -v pop="$pop" '$2==pop {print $1}' "$PANEL_FILE" > "$SAMPLES_DIR/all_${pop}.txt"
    n=$(wc -l < "$SAMPLES_DIR/all_${pop}.txt")
    echo "  $pop: $n muestras disponibles"
done
echo ""

select_samples() {
    shuf -n "$2" --random-source=<(yes $3) "$1" > "$4"
}

process_chromosome() {
    local chr=$1 run=$2 run_dir=$3 output_prefix=$4
    local vcf="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    for pop in "${POPULATIONS[@]}"; do
        samples="$run_dir/samples_${pop}.txt"
        bcftools view -S "$samples" -Oz -o "${output_prefix}_${pop}.vcf.gz" "$vcf" 2>/dev/null
        bcftools +fill-tags "${output_prefix}_${pop}.vcf.gz" -Oz -o "${output_prefix}_${pop}_AF.vcf.gz" -- -t AF 2>/dev/null
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
    
    rm -f "${output_prefix}"_{IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM}{,_AF}.vcf.gz \
          "${output_prefix}"_{IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM}.freq
    
    echo "    [✓] Chr$chr"
}

# BUCLE PRINCIPAL
for run in {1..5}; do
    echo ""
    echo "========== CORRIDA $run/$N_RUNS =========="
    echo "Inicio: $(date)"
    
    run_dir="$WORK_DIR/panel_N11_runs/run${run}"
    mkdir -p "$run_dir"
    
    echo "Seleccionando muestras..."
    for pop in "${POPULATIONS[@]}"; do
        select_samples "$SAMPLES_DIR/all_${pop}.txt" $N_SAMPLES $run "$run_dir/samples_${pop}.txt"
    done
    
    echo "Procesando cromosomas (${PARALLEL_CHR} paralelos)..."
    pids=()
    for chr in {1..22}; do
        process_chromosome $chr $run "$run_dir" "$run_dir/chr${chr}" &
        pids+=($!)
        
        if [[ ${#pids[@]} -ge $PARALLEL_CHR ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
    for pid in "${pids[@]}"; do wait "$pid"; done
    
    echo "Concatenando..."
    {
        echo -e "#chrom\tposition\trsid\tA1\tA2\tIBS\tCEU\tYRI\tESN\tGWD\tLWK\tPEL\tMXL\tPUR\tFIN\tCLM"
        for chr in {1..22}; do
            tail -n +2 "$run_dir/chr${chr}.freq" 2>/dev/null || true
        done
    } > "$OUTPUT_DIR/PANEL_N11_run${run}.txt"
    
    # Filtrar MAF > 0.1 en al menos una población
    echo "Aplicando filtro MAF > 0.1..."
    awk 'BEGIN{FS=OFS="\t"} 
         NR==1 {print; next}
         {
            has_informative = 0
            for(i=6; i<=16; i++) {
                if($i > 0.1 && $i < 0.9) has_informative = 1
            }
            if(has_informative) print
         }' "$OUTPUT_DIR/PANEL_N11_run${run}.txt" > "$OUTPUT_DIR/PANEL_N11_run${run}_extreme.txt"
    
    rm -rf "$run_dir/chr"*.freq
    
    nsnps=$(awk 'NR>1' "$OUTPUT_DIR/PANEL_N11_run${run}_extreme.txt" | wc -l)
    size=$(du -h "$OUTPUT_DIR/PANEL_N11_run${run}_extreme.txt" | cut -f1)
    echo "✓ Run$run: $nsnps SNPs ($size) - $(date)"
done

echo ""
echo "=============================================="
echo "✓ 5 PANELES N11 FILTRADOS COMPLETADOS"
echo "=============================================="
ls -lh "$OUTPUT_DIR"/PANEL_N11_run*_extreme.txt
