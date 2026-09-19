#!/bin/bash
set -euo pipefail

RESULTS_DIR="$HOME/Documentos/Ancestria/results/genetic_distances"
FREQ_DIR="$HOME/Documentos/Ancestria/work/freq_genome"
CALC="$HOME/Documentos/Ancestria/work/calc_distances.py"

echo "Calculando distancias N11..."

# Las frecuencias N11 ya fueron extraídas
# Solo ejecutar cálculos

pids=()
for sample in POOL1 POOL2 HOSPITAL; do
    for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n11.tsv"
        ref="$FREQ_DIR/${pop}_n11.tsv"
        
        # Verificar ref existe
        if [[ ! -s "$ref" ]]; then
            echo "  [!] Extrayendo $pop..."
            awk -v p="$pop" 'BEGIN{FS=OFS="\t"} 
                NR==1{for(i=1;i<=NF;i++)if($i==p)c=i;next} 
                c{print $1,$2,$4,$5,$c}' \
                "$HOME/Documentos/Ancestria/refpanels/FINAL_FREQUENCIESN11.txt" > "$ref"
        fi
        
        python3 "$CALC" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$ref" \
            --sample "$sample" \
            --pop "$pop" \
            --out "$out" 2>&1 &
        
        pids+=($!)
        
        # Limitar a 6 procesos simultáneos
        if [[ ${#pids[@]} -ge 6 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

for pid in "${pids[@]}"; do wait "$pid"; done

# Consolidar
summary="$RESULTS_DIR/summary_distances_n11.tsv"
echo -e "Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya" > "$summary"
for f in "$RESULTS_DIR"/*_vs_*_n11.tsv; do
    [[ -s "$f" ]] && tail -1 "$f" >> "$summary"
done

echo ""
echo "=== N11 ==="
cat "$summary" | column -t

echo ""
echo "✅ Completado"
