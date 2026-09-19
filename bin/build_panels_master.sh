#!/bin/bash
# ==============================================================================
# CONSTRUCCIÓN DE PANELES: 10 CORRIDAS
# ==============================================================================
# Panel A: 3 Superpoblaciones (AFR, EUR, AMR) × 10 corridas
# Panel B: N11 (11 poblaciones específicas) × 10 corridas
# ==============================================================================

set -euo pipefail

# Configuración
PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/data/vcf"
PANEL_FILE="$PROJ_DIR/data/panel/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR="$PROJ_DIR/refpanels"
WORK_DIR="$PROJ_DIR/work"

source "$PROJ_DIR/env.sh"

echo "=============================================="
echo "CONSTRUCCIÓN DE PANELES - 10 CORRIDAS"
echo "=============================================="
echo ""
echo "Verificando archivos necesarios..."

# Verificar que existan los VCFs
n_vcfs=$(ls "$VCF_DIR"/ALL.chr*.vcf.gz 2>/dev/null | wc -l)

if [ "$n_vcfs" -lt 22 ]; then
    echo "ERROR: Faltan VCFs del genoma completo"
    echo "Encontrados: $n_vcfs/22"
    echo ""
    echo "Ejecuta primero: ./download_1kgp_phase3.sh"
    exit 1
fi

echo "✓ VCFs encontrados: $n_vcfs/22"

# Verificar panel file
if [ ! -f "$PANEL_FILE" ]; then
    echo "ERROR: No se encuentra $PANEL_FILE"
    exit 1
fi

echo "✓ Panel file encontrado"
echo ""

# Contar muestras disponibles
n_afr=$(awk '$3=="AFR"' "$PANEL_FILE" | wc -l)
n_eur=$(awk '$3=="EUR"' "$PANEL_FILE" | wc -l)
n_amr=$(awk '$3=="AMR"' "$PANEL_FILE" | wc -l)

echo "Muestras disponibles:"
echo "  AFR: $n_afr"
echo "  EUR: $n_eur"
echo "  AMR: $n_amr"
echo ""

# Verificar suficientes muestras
if [ "$n_afr" -lt 100 ] || [ "$n_eur" -lt 100 ] || [ "$n_amr" -lt 100 ]; then
    echo "ERROR: Necesitas al menos 100 muestras de cada superpoblación"
    exit 1
fi

echo "✓ Suficientes muestras para 10 corridas"
echo ""

# Menú de selección
echo "¿Qué paneles quieres construir?"
echo ""
echo "1) Solo Panel A (3 superpoblaciones) - ~6-8 horas"
echo "2) Solo Panel B (N11, 11 poblaciones) - ~8-10 horas"  
echo "3) Ambos paneles (recomendado) - ~14-18 horas"
echo ""
read -p "Selecciona opción (1/2/3): " option

case $option in
    1)
        echo ""
        echo "Construyendo Panel A (3 superpoblaciones)..."
        bash "$PROJ_DIR/bin/build_panel_3way_10runs.sh"
        ;;
    2)
        echo ""
        echo "Construyendo Panel B (N11)..."
        bash "$PROJ_DIR/bin/build_panel_N11_10runs.sh"
        ;;
    3)
        echo ""
        echo "Construyendo ambos paneles..."
        bash "$PROJ_DIR/bin/build_panel_3way_10runs.sh"
        bash "$PROJ_DIR/bin/build_panel_N11_10runs.sh"
        ;;
    *)
        echo "Opción inválida"
        exit 1
        ;;
esac

echo ""
echo "=============================================="
echo "CONSTRUCCIÓN COMPLETADA"
echo "=============================================="
echo ""
echo "Paneles creados en: $OUTPUT_DIR"
echo ""
echo "Siguiente paso: Ejecutar iAdmix con los paneles"
echo "  ./run_iadmix_all.sh"
