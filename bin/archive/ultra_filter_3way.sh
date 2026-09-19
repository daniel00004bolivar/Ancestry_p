#!/bin/bash
for run in {1..10}; do
    input="$HOME/Documentos/Ancestria/refpanels/PANEL_3WAY_run${run}_final.txt"
    output="$HOME/Documentos/Ancestria/refpanels/PANEL_3WAY_run${run}_ultra.txt"
    
    awk 'BEGIN{FS=OFS="\t"} 
         NR==1 {print; next}
         ($6>0.05 || $7>0.05 || $8>0.05) && ($6<0.95 || $7<0.95 || $8<0.95)
        ' "$input" > "$output"
    
    echo "Run $run: $(wc -l < $output) SNPs"
done
