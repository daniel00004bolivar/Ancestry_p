#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
RESULTS_3WAY="$PROJ_DIR/results/iadmix_10runs_3way"
RESULTS_N11="$PROJ_DIR/results/iadmix_5runs_N11_FULL"
OUTPUT_DIR="$PROJ_DIR/results/distances"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "CALCULANDO DISTANCIAS GENÉTICAS"
echo "=========================================="
echo ""

# ============================================
# Panel 3-WAY
# ============================================

echo "=== PANEL 3-WAY ===" > "$OUTPUT_DIR/distances_3way.txt"
echo "" >> "$OUTPUT_DIR/distances_3way.txt"

RUNS_3WAY=(1 2 3 6 7)

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:" >> "$OUTPUT_DIR/distances_3way.txt"
    
    distances=()
    for ((i=0; i<${#RUNS_3WAY[@]}; i++)); do
        run1=${RUNS_3WAY[$i]}
        for ((j=i+1; j<${#RUNS_3WAY[@]}; j++)); do
            run2=${RUNS_3WAY[$j]}
            
            out1="$RESULTS_3WAY/${sample}_run${run1}.ancestry.out"
            out2="$RESULTS_3WAY/${sample}_run${run2}.ancestry.out"
            
            if [[ -s "$out1" && -s "$out2" ]]; then
                # Extraer AFR, EUR, AMR
                v1_afr=$(grep "FINAL_NZ_PROPS" "$out1" | grep -oP "AFR:[0-9.]+" | cut -d: -f2)
                v1_eur=$(grep "FINAL_NZ_PROPS" "$out1" | grep -oP "EUR:[0-9.]+" | cut -d: -f2)
                v1_amr=$(grep "FINAL_NZ_PROPS" "$out1" | grep -oP "AMR:[0-9.]+" | cut -d: -f2)
                
                v2_afr=$(grep "FINAL_NZ_PROPS" "$out2" | grep -oP "AFR:[0-9.]+" | cut -d: -f2)
                v2_eur=$(grep "FINAL_NZ_PROPS" "$out2" | grep -oP "EUR:[0-9.]+" | cut -d: -f2)
                v2_amr=$(grep "FINAL_NZ_PROPS" "$out2" | grep -oP "AMR:[0-9.]+" | cut -d: -f2)
                
                # Distancia euclidea
                dist=$(awk -v a1="$v1_afr" -v e1="$v1_eur" -v m1="$v1_amr" \
                           -v a2="$v2_afr" -v e2="$v2_eur" -v m2="$v2_amr" \
                           'BEGIN{printf "%.6f", sqrt((a1-a2)^2 + (e1-e2)^2 + (m1-m2)^2)}')
                
                echo "  Run$run1 vs Run$run2: $dist" >> "$OUTPUT_DIR/distances_3way.txt"
                distances+=("$dist")
            fi
        done
    done
    
    # Promedio
    if [[ ${#distances[@]} -gt 0 ]]; then
        avg=$(printf '%s\n' "${distances[@]}" | awk '{sum+=$1; sum2+=$1*$1; n++} END{printf "%.6f ± %.6f", sum/n, sqrt((sum2-n*(sum/n)^2)/(n-1))}')
        echo "  Promedio: $avg" >> "$OUTPUT_DIR/distances_3way.txt"
    fi
    echo "" >> "$OUTPUT_DIR/distances_3way.txt"
done

# ============================================
# Panel N11
# ============================================

echo "=== PANEL N11 ===" > "$OUTPUT_DIR/distances_N11.txt"
echo "" >> "$OUTPUT_DIR/distances_N11.txt"

POPS_N11=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:" >> "$OUTPUT_DIR/distances_N11.txt"
    
    distances=()
    for ((i=1; i<=4; i++)); do
        for ((j=i+1; j<=5; j++)); do
            out1="$RESULTS_N11/${sample}_run${i}.ancestry.out"
            out2="$RESULTS_N11/${sample}_run${j}.ancestry.out"
            
            if [[ -s "$out1" && -s "$out2" ]]; then
                # Extraer valores
                vals1=""
                vals2=""
                for pop in "${POPS_N11[@]}"; do
                    v1=$(grep "FINAL_NZ_PROPS" "$out1" | grep -oP "${pop}:[0-9.]+" | cut -d: -f2 || echo "0")
                    v2=$(grep "FINAL_NZ_PROPS" "$out2" | grep -oP "${pop}:[0-9.]+" | cut -d: -f2 || echo "0")
                    vals1="$vals1 $v1"
                    vals2="$vals2 $v2"
                done
                
                # Distancia euclidea
                dist=$(echo "$vals1" "$vals2" | awk '{
                    sum=0
                    n=NF/2
                    for(i=1; i<=n; i++) {
                        sum += ($(i) - $(i+n))^2
                    }
                    printf "%.6f", sqrt(sum)
                }')
                
                echo "  Run$i vs Run$j: $dist" >> "$OUTPUT_DIR/distances_N11.txt"
                distances+=("$dist")
            fi
        done
    done
    
    # Promedio
    if [[ ${#distances[@]} -gt 0 ]]; then
        avg=$(printf '%s\n' "${distances[@]}" | awk '{sum+=$1; sum2+=$1*$1; n++} END{if(n>1) printf "%.6f ± %.6f", sum/n, sqrt((sum2-n*(sum/n)^2)/(n-1)); else printf "%.6f", sum/n}')
        echo "  Promedio: $avg" >> "$OUTPUT_DIR/distances_N11.txt"
    fi
    echo "" >> "$OUTPUT_DIR/distances_N11.txt"
done

echo ""
echo "=========================================="
echo "✓ DISTANCIAS CALCULADAS"
echo "=========================================="
echo ""
echo "Panel 3-WAY:"
cat "$OUTPUT_DIR/distances_3way.txt"
echo ""
echo "Panel N11:"
cat "$OUTPUT_DIR/distances_N11.txt"
