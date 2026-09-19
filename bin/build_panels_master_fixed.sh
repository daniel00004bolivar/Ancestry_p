#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
BIN_DIR="$PROJ_DIR/bin"

echo "¿Qué paneles quieres construir?"
echo "1) Solo Panel A (3 superpoblaciones)"
echo "2) Solo Panel B (N11)"  
echo "3) Ambos paneles"
read -p "Opción (1/2/3): " option

case $option in
    1) bash "$BIN_DIR/build_panel_3way_10runs.sh" ;;
    2) bash "$BIN_DIR/build_panel_N11_10runs.sh" ;;
    3) 
        bash "$BIN_DIR/build_panel_3way_10runs.sh"
        bash "$BIN_DIR/build_panel_N11_10runs.sh"
        ;;
    *) echo "Opción inválida"; exit 1 ;;
esac
