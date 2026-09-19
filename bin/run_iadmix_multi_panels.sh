#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
RESULTS_DIR="$PROJ_DIR/results/iadmix_multi_panels"

mkdir -p "$RESULTS_DIR"

# Paneles a probar (DEBEN EXISTIR)
declare -A PANELS=(
    [N3]="$PROJ_DIR/refpanels/PANEL_3WAY_run1_extreme.txt"
    [N5]="$PROJ_DIR/refpanels/FINAL_FREQUENCIES_GNOMAD3.txt"
    [N11]="$PROJ_DIR/refpanels/FINAL_FREQUENCIESN11.txt"
)

# BAMs
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
echo "iADMIX MULTI-PANEL: N3, N5, N11"
echo "3 paneles × 3 muestras = 9 análisis"
echo "Paralelizado: 3 muestras simultáneas"
echo "=========================================="
echo ""

# Verificar paneles
for panel_name in "${!PANELS[@]}"; do
    panel="${PANELS[$panel_name]}"
    if [[ ! -f "$panel" ]]; then
        echo "ERROR: Panel $panel_name no existe: $panel"
        exit 1
    fi
    nsnps=$(awk 'NR>1' "$panel" | wc -l)
    echo "  ✓ $panel_name: $nsnps SNPs"
done
echo ""

# Ejecutar análisis
for panel_name in N3 N5 N11; do
    panel="${PANELS[$panel_name]}"
    
    echo "=========================================="
    echo "Panel: $panel_name"
    echo "=========================================="
    
    # Correr las 3 muestras EN PARALELO
    pids=()
    for sample in POOL1 POOL2 HOSPITAL; do
        bam="${BAMS[$sample]}"
        poolsize="${POOLSIZES[$sample]}"
        output="$RESULTS_DIR/${sample}_${panel_name}"
        
        if [[ -s "${output}.ancestry.out" ]]; then
            echo "  [=] $sample ya existe"
            continue
        fi
        
        echo "  [..] $sample (n=$poolsize) - inicio: $(date +%H:%M:%S)"
        
        (
            python2 "$IADMIX/runancestry.py" \
                -f "$panel" \
                --bam "$bam" \
                -p "$poolsize" \
                -o "$output" \
                --path "$IADMIX" \
                > "${output}.log" 2>&1
            
            if [[ -s "${output}.ancestry.out" ]]; then
                echo "  [✓] $sample completado: $(date +%H:%M:%S)"
            else
                echo "  [✗] $sample FALLÓ"
            fi
        ) &
        pids+=($!)
    done
    
    # Esperar a que terminen las 3
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
    
    echo ""
done

# Resumen
echo "=========================================="
echo "RESUMEN DE RESULTADOS"
echo "=========================================="
echo ""

for panel_name in N3 N5 N11; do
    echo "Panel: $panel_name"
    for sample in POOL1 POOL2 HOSPITAL; do
        out="$RESULTS_DIR/${sample}_${panel_name}.ancestry.out"
        if [[ -s "$out" ]]; then
            result=$(grep "ADMIX_PROP\|FINAL_NZ_PROPS" "$out" | tail -1)
            echo "  $sample: $result"
        else
            echo "  $sample: FALLO"
        fi
    done
    echo ""
done
