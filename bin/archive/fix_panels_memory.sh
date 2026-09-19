#!/bin/bash
set -euo pipefail

REFPANELS="$HOME/Documentos/Ancestria/refpanels"

echo "=============================================="
echo "REDUCIR PANELES PARA EVITAR OOM"
echo "=============================================="

for run in {1..10}; do
    input="$REFPANELS/PANEL_3WAY_run${run}_clean.txt"
    output="$REFPANELS/PANEL_3WAY_run${run}_final.txt"
    
    if [ ! -f "$input" ]; then
        echo "Run $run: no existe, saltando"
        continue
    fi
    
    echo "Procesando run $run..."
    
    # Filtrar:
    # 1. Solo SNPs bialélicos (ya hecho en _clean)
    # 2. MAF > 0.01 en AL MENOS una población
    # 3. NO monomórficos (frecuencia 0 o 1 en TODAS las poblaciones)
    
    awk 'BEGIN{FS=OFS="\t"}
         NR==1 {print; next}
         {
           afr=$6; eur=$7; amr=$8
           
           # Eliminar monomórficos (0 o 1 en todas)
           if ((afr==0 || afr==1) && (eur==0 || eur==1) && (amr==0 || amr==1)) next
           
           # Mantener si MAF > 0.01 en alguna población
           if (afr>0.01 && afr<0.99) {print; next}
           if (eur>0.01 && eur<0.99) {print; next}
           if (amr>0.01 && amr<0.99) {print; next}
         }
        ' "$input" > "$output"
    
    orig=$(wc -l < "$input")
    filt=$(wc -l < "$output")
    diff=$((orig - filt))
    pct=$(awk "BEGIN{printf \"%.1f\", ($diff/$orig)*100}")
    
    echo "  Original: $orig SNPs"
    echo "  Filtrado: $filt SNPs"
    echo "  Removido: $diff SNPs ($pct%)"
    echo ""
done

echo "=============================================="
echo "✓ PANELES OPTIMIZADOS"
echo "=============================================="
echo ""
ls -lh "$REFPANELS"/PANEL_3WAY_run*_final.txt
