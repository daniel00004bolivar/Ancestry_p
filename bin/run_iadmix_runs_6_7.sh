#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_10runs_3way"

declare -A BAMS=(
    [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
)

declare -A POOLSIZES=([POOL1]=50 [POOL2]=50 [HOSPITAL]=200)

echo "=========================================="
echo "iADMIX 3-WAY: Runs 6 y 7"
echo "=========================================="

for run in 6 7; do
    panel="$PANELS_DIR/PANEL_3WAY_run${run}_extreme.txt"
    
    echo ""
    echo "=========================================="
    echo "CORRIDA $run"
    echo "Inicio: $(date)"
    echo "=========================================="
    
    # Correr 3 muestras en paralelo
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
            else
                echo "  [✗] $sample run$run FALLÓ"
            fi
        ) &
        pids+=($!)
    done
    
    for pid in "${pids[@]}"; do wait "$pid"; done
    echo "  Corrida $run completada: $(date)"
done

echo ""
echo "=========================================="
echo "✓ Runs 6 y 7 completados"
echo "=========================================="
