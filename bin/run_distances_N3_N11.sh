#!/bin/bash
set -euo pipefail

RESULTS_DIR="$HOME/Documentos/Ancestria/results/genetic_distances"
FREQ_DIR="$HOME/Documentos/Ancestria/work/freq_genome"
CALC_SCRIPT="$HOME/Documentos/Ancestria/work/calc_distances.py"

echo "=============================================="
echo "Calculando distancias N3 y N11"
echo "=============================================="

pids=()

# N3
echo ""
echo "Panel N3..."
for sample in POOL1 POOL2 HOSPITAL; do
    for pop in AFR EUR AMR; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n3.tsv"
        
        if [[ -s "$out" ]]; then
            echo "  [=] $sample vs $pop (N3)"
            continue
        fi
        
        python3 "$CALC_SCRIPT" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$FREQ_DIR/${pop}_n3.tsv" \
            --sample "$sample" \
            --pop "$pop" \
            --out "$out" 2>&1 | grep -v "Cargando" &
        
        pids+=($!)
        
        if [[ ${#pids[@]} -ge 6 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

# N11
echo ""
echo "Panel N11..."
for sample in POOL1 POOL2 HOSPITAL; do
    for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n11.tsv"
        
        if [[ -s "$out" ]]; then
            echo "  [=] $sample vs $pop (N11)"
            continue
        fi
        
        # Verificar que existe el archivo de frecuencias
        ref_file="$FREQ_DIR/${pop}_n11.tsv"
        if [[ ! -s "$ref_file" ]]; then
            echo "  [!] Falta $ref_file - extrayendo..."
            awk -v p="$pop" 'BEGIN{FS=OFS="\t"} 
                NR==1{for(i=1;i<=NF;i++)if($i==p)c=i;next} 
                c{print $1,$2,$4,$5,$c}' \
                "$HOME/Documentos/Ancestria/refpanels/FINAL_FREQUENCIESN11.txt" > "$ref_file"
        fi
        
        python3 "$CALC_SCRIPT" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$ref_file" \
            --sample "$sample" \
            --pop "$pop" \
            --out "$out" 2>&1 | grep -v "Cargando" &
        
        pids+=($!)
        
        if [[ ${#pids[@]} -ge 6 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

for pid in "${pids[@]}"; do wait "$pid"; done

# Consolidar
echo ""
echo "=============================================="
echo "Consolidando resultados..."
echo "=============================================="

for panel in n3 n11; do
    summary="$RESULTS_DIR/summary_distances_${panel}.tsv"
    echo -e "Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya" > "$summary"
    
    for f in "$RESULTS_DIR"/*_vs_*_${panel}.tsv; do
        [[ -s "$f" ]] && tail -1 "$f" >> "$summary"
    done
    
    echo ""
    echo "=== ${panel^^} ==="
    cat "$summary" | column -t
done

echo ""
echo "✅ Completado"
ls -lh "$RESULTS_DIR"/summary_distances_*.tsv
