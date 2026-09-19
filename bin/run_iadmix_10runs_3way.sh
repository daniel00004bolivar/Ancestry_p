#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_10runs_3way"

mkdir -p "$RESULTS_DIR"

declare -A BAMS=(
    [POOL1]="$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    [POOL2]="$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    [HOSPITAL]="$PROJ_DIR/data/bam/HOSPITAL.bam"
)

declare -A POOLSIZES=(
    [POOL1]=50
    [POOL2]=50
    [HOSPITAL]=200
)

echo "=========================================="
echo "iADMIX 3-WAY: 10 CORRIDAS"
echo "Objetivo: Medir robustez (media ± SD)"
echo "10 paneles × 3 muestras = 30 análisis"
echo "=========================================="
echo ""

# Verificar paneles
for run in {1..10}; do
    panel="$PANELS_DIR/PANEL_3WAY_run${run}_extreme.txt"
    if [[ ! -f "$panel" ]]; then
        echo "ERROR: Falta panel run${run}"
        exit 1
    fi
done
echo "✓ 10 paneles verificados"
echo ""

# Ejecutar 10 corridas
for run in {1..10}; do
    panel="$PANELS_DIR/PANEL_3WAY_run${run}_extreme.txt"
    
    echo "=========================================="
    echo "CORRIDA $run/10"
    echo "Panel: PANEL_3WAY_run${run}_extreme.txt"
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
            
            echo "  [✓] $sample run$run: $(date +%H:%M:%S)"
        ) &
        pids+=($!)
    done
    
    for pid in "${pids[@]}"; do wait "$pid"; done
    echo "  Corrida $run completada: $(date)"
    echo ""
done

# Calcular estadísticas
echo "=========================================="
echo "CALCULANDO MEDIA ± SD"
echo "=========================================="

for sample in POOL1 POOL2 HOSPITAL; do
    echo ""
    echo "$sample:"
    
    for pop in AFR EUR AMR; do
        vals=()
        for run in {1..10}; do
            out="$RESULTS_DIR/${sample}_run${run}.ancestry.out"
            if [[ -s "$out" ]]; then
                val=$(grep "FINAL_NZ_PROPS\|ADMIX_PROP" "$out" | tail -1 | grep -oP "${pop}:[0-9.]+" | cut -d: -f2)
                [[ -n "$val" ]] && vals+=("$val")
            fi
        done
        
        if [[ ${#vals[@]} -gt 0 ]]; then
            stats=$(printf '%s\n' "${vals[@]}" | awk '{sum+=$1; sum2+=$1*$1; n++} END{mean=sum/n; sd=sqrt((sum2-n*mean*mean)/(n-1)); printf "%.4f ± %.4f", mean, sd}')
            echo "  $pop: $stats"
        fi
    done
done

echo ""
echo "=========================================="
echo "✓ COMPLETADO"
echo "Resultados: $RESULTS_DIR"
echo "=========================================="
