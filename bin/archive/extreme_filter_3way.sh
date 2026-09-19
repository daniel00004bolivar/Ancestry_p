#!/bin/bash
# MAF > 0.1 para ~3M SNPs
for run in {1..10}; do
    input="$HOME/Documentos/Ancestria/refpanels/PANEL_3WAY_run${run}_final.txt"
    output="$HOME/Documentos/Ancestria/refpanels/PANEL_3WAY_run${run}_extreme.txt"
    
    awk 'BEGIN{FS=OFS="\t"} 
         NR==1 {print; next}
         ($6>0.1 || $7>0.1 || $8>0.1) && ($6<0.9 || $7<0.9 || $8<0.9)
        ' "$input" > "$output"
    
    orig=$(wc -l < "$input")
    filt=$(wc -l < "$output")
    echo "Run $run: $orig → $filt SNPs"
done

echo ""
echo "Tamaño archivos:"
ls -lh "$HOME/Documentos/Ancestria/refpanels"/PANEL_3WAY_run*_extreme.txt
