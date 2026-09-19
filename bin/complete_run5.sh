#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_10runs_3way"

panel="$PANELS_DIR/PANEL_3WAY_run5_extreme.txt"
sample="POOL1"
bam="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
poolsize=50
output="$RESULTS_DIR/${sample}_run5"

echo "=========================================="
echo "Completando POOL1 run5"
echo "=========================================="

python2 "$IADMIX/runancestry.py" \
    -f "$panel" \
    --bam "$bam" \
    -p "$poolsize" \
    -o "$output" \
    --path "$IADMIX" \
    2>&1 | tee "${output}.log"

if [[ -s "${output}.ancestry.out" ]]; then
    echo ""
    echo "✓ POOL1 run5 completado"
    grep "FINAL_NZ_PROPS\|ADMIX_PROP" "${output}.ancestry.out" | tail -1
else
    echo "✗ Error - revisar ${output}.log"
fi
