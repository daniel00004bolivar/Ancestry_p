#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_3way_10runs"

mkdir -p "$RESULTS_DIR"

POOLS=("POOL1" "POOL2" "HOSPITAL")
BAMS=(
    "$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    "$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    "$PROJ_DIR/data/bam/HOSPITAL.bam"
)
POOLSIZES=(50 50 200)  # ← CORRECCIÓN

echo "=========================================="
echo "iADMIX: 30 ANÁLISIS (10 paneles × 3 pools)"
echo "=========================================="
echo ""

total=0
for run in {1..10}; do
    panel="$PANELS_DIR/PANEL_3WAY_run${run}_ultra.txt"
    
    if [ ! -f "$panel" ]; then
        echo "ERROR: No existe $panel"
        exit 1
    fi
    
    for i in {0..2}; do
        pool="${POOLS[$i]}"
        bam="${BAMS[$i]}"
        poolsize="${POOLSIZES[$i]}"
        output="$RESULTS_DIR/${pool}_3WAY_run${run}"
        
        total=$((total + 1))
        echo "[$total/30] Panel run$run - $pool (n=$poolsize)"
        
        python2 "$IADMIX/runancestry.py" \
            -f "$panel" \
            --bam "$bam" \
            -p "$poolsize" \
            -o "$output" \
            --path "$IADMIX"
        
        echo "  ✓ Completado"
        echo ""
    done
done

echo ""
echo "=========================================="
echo "✓ 30 ANÁLISIS COMPLETADOS"
echo "=========================================="
