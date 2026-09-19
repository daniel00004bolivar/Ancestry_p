#!/bin/bash
# Limpia paneles: solo SNPs bialélicos (A,T,C,G)
set -euo pipefail

REFPANELS="$HOME/Documentos/Ancestria/refpanels"

echo "=============================================="
echo "LIMPIEZA: Solo SNPs bialélicos"
echo "=============================================="
echo ""

for run in {1..10}; do
    panel="$REFPANELS/PANEL_3WAY_run${run}.txt"
    clean="$REFPANELS/PANEL_3WAY_run${run}_clean.txt"
    
    if [ ! -f "$panel" ]; then
        echo "⏳ Run $run: pendiente"
        continue
    fi
    
    echo "Procesando Run $run..."
    
    # Filtrar: solo A,T,C,G bialélicos
    awk 'BEGIN{FS=OFS="\t"}
         NR==1 {print; next}
         $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ && $4 != $5
        ' "$panel" > "$clean"
    
    orig=$(wc -l < "$panel")
    filt=$(wc -l < "$clean")
    diff=$((orig - filt))
    pct=$(awk "BEGIN{printf \"%.1f\", ($diff/$orig)*100}")
    
    echo "  Original: $orig SNPs"
    echo "  Filtrado: $filt SNPs"
    echo "  Removido: $diff SNPs ($pct%)"
    echo ""
done

echo "=============================================="
echo "✓ LIMPIEZA COMPLETADA"
echo "=============================================="
echo ""
echo "Paneles limpios en:"
ls -lh "$REFPANELS"/PANEL_3WAY_run*_clean.txt
