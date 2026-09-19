#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_5runs_N11_FULL"

mkdir -p "$RESULTS_DIR"

declare -A BAMS=(
    [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
)

declare -A POOLSIZES=([POOL1]=50 [POOL2]=50 [HOSPITAL]=200)

echo "=========================================="
echo "iADMIX N11 COMPLETO: 5 CORRIDAS"
echo "81M SNPs × 5 paneles × 3 muestras"
echo "=========================================="

for run in {1..5}; do
    panel="$PANELS_DIR/PANEL_N11_FULL_run${run}.txt"
    
    echo ""
    echo "=========================================="
    echo "CORRIDA $run/5"
    echo "Inicio: $(date)"
    echo "=========================================="
    
    pids=()
    for sample in POOL1 POOL2 HOSPITAL; do
        bam="${BAMS[$sample]}"
        poolsize="${POOLSIZES[$sample]}"
        output="$RESULTS_DIR/${sample}_run${run}"
        
        if [[ -s "${output}.ancestry.out" ]]; then
            echo "  [=] $sample run$run ya existe"
            continue
        fi
        
        echo "  [..] $sample (n=$poolsize)"
        
        (
            python2 "$IADMIX/runancestry.py" \
                -f "$panel" \
                --bam "$bam" \
                -p "$poolsize" \
                -o "$output" \
                --path "$IADMIX" \
                > "${output}.log" 2>&1
            
            if [[ -s "${output}.ancestry.out" ]]; then
                echo "  [✓] $sample run$run: $(date +%H:%M:%S)"
            fi
        ) &
        pids+=($!)
    done
    
    for pid in "${pids[@]}"; do wait "$pid"; done
    echo "  Corrida $run completada: $(date)"
done

echo ""
echo "=========================================="
echo "CALCULANDO ESTADÍSTICAS N11..."
echo "=========================================="

POPS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

for sample in POOL1 POOL2 HOSPITAL; do
    echo ""
    echo "$sample:"
    for pop in "${POPS[@]}"; do
        vals=()
        for run in {1..5}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            [[ -s "$out" ]] && vals+=($(grep "FINAL_NZ_PROPS" "$out" | grep -oP "${pop}:[0-9.]+" | cut -d: -f2 || echo ""))
        done
        
        [[ ${#vals[@]} -eq 5 ]] && printf "  %-4s: %s\n" "$pop" "$(printf '%s\n' "${vals[@]}" | awk '{sum+=$1; sum2+=$1*$1; n++; v[n]=$1} END{printf "%.4f ± %.4f  [%.4f, %.4f, %.4f, %.4f, %.4f]", sum/n, sqrt((sum2-n*(sum/n)^2)/(n-1)), v[1], v[2], v[3], v[4], v[5]}')"
    done
done

echo ""
echo "✓ COMPLETADO"
