#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_5runs_N11"

mkdir -p "$RESULTS_DIR"

declare -A BAMS=(
    [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
)

declare -A POOLSIZES=([POOL1]=50 [POOL2]=50 [HOSPITAL]=200)

echo "=========================================="
echo "iADMIX N11: 5 CORRIDAS"
echo "11 poblaciones × 5 corridas × 3 muestras"
echo "=========================================="

for run in {1..5}; do
    panel="$PANELS_DIR/PANEL_N11_run${run}_extreme.txt"
    
    echo ""
    echo "=========================================="
    echo "CORRIDA $run/5"
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
echo "CALCULANDO ESTADÍSTICAS N11..."
echo "=========================================="
echo ""

POPS=("IBS" "CEU" "YRI" "ESN" "GWD" "LWK" "PEL" "MXL" "PUR" "FIN" "CLM")

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:"
    for pop in "${POPS[@]}"; do
        echo -n "  $pop: "
        vals=()
        for run in {1..5}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            if [[ -s "$out" ]]; then
                val=$(grep "FINAL_NZ_PROPS" "$out" | grep -oP "${pop}:[0-9.]+" | cut -d: -f2)
                [[ -n "$val" ]] && vals+=("$val")
            fi
        done
        
        if [[ ${#vals[@]} -eq 5 ]]; then
            printf '%s\n' "${vals[@]}" | awk '{
                sum+=$1; sum2+=$1*$1; n++
                v[n]=$1
            } END {
                mean=sum/n
                sd=sqrt((sum2-n*mean*mean)/(n-1))
                printf "%.4f ± %.4f\n", mean, sd
            }'
        fi
    done
    echo ""
done

echo ""
echo "=========================================="
echo "Agregación 3-WAY desde N11:"
echo "=========================================="

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:"
    
    # AFR = YRI + ESN + GWD + LWK
    echo -n "  AFR (YRI+ESN+GWD+LWK): "
    for comp in YRI ESN GWD LWK; do
        vals=()
        for run in {1..5}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            [[ -s "$out" ]] && vals+=($(grep "FINAL_NZ_PROPS" "$out" | grep -oP "${comp}:[0-9.]+" | cut -d: -f2))
        done
        [[ ${#vals[@]} -eq 5 ]] && printf '%s\n' "${vals[@]}"
    done | awk '{for(i=1;i<=5;i++) sum[i]+=$1} END{for(i=1;i<=5;i++) v[i]=sum[i]; n=5; for(i in v){s+=v[i]; s2+=v[i]*v[i]} printf "%.4f ± %.4f\n", s/n, sqrt((s2-n*(s/n)^2)/(n-1))}'
    
    # EUR = IBS + CEU + FIN
    echo -n "  EUR (IBS+CEU+FIN): "
    for comp in IBS CEU FIN; do
        vals=()
        for run in {1..5}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            [[ -s "$out" ]] && vals+=($(grep "FINAL_NZ_PROPS" "$out" | grep -oP "${comp}:[0-9.]+" | cut -d: -f2))
        done
        [[ ${#vals[@]} -eq 5 ]] && printf '%s\n' "${vals[@]}"
    done | awk '{for(i=1;i<=5;i++) sum[i]+=$1} END{for(i=1;i<=5;i++) v[i]=sum[i]; n=5; for(i in v){s+=v[i]; s2+=v[i]*v[i]} printf "%.4f ± %.4f\n", s/n, sqrt((s2-n*(s/n)^2)/(n-1))}'
    
    # AMR = PEL + MXL + PUR + CLM
    echo -n "  AMR (PEL+MXL+PUR+CLM): "
    for comp in PEL MXL PUR CLM; do
        vals=()
        for run in {1..5}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            [[ -s "$out" ]] && vals+=($(grep "FINAL_NZ_PROPS" "$out" | grep -oP "${comp}:[0-9.]+" | cut -d: -f2))
        done
        [[ ${#vals[@]} -eq 5 ]] && printf '%s\n' "${vals[@]}"
    done | awk '{for(i=1;i<=5;i++) sum[i]+=$1} END{for(i=1;i<=5;i++) v[i]=sum[i]; n=5; for(i in v){s+=v[i]; s2+=v[i]*v[i]} printf "%.4f ± %.4f\n", s/n, sqrt((s2-n*(s/n)^2)/(n-1))}'
    
    echo ""
done

echo "✓ COMPLETADO"
