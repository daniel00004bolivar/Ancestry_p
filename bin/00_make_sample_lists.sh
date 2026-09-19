#!/usr/bin/env bash
set -euo pipefail

# Uso:
#   bash bin/00_make_sample_lists.sh <POPS_FILE>
#
# Lee POPS_FILE (p.ej. config/pops_n5.txt) y crea data/samples/<POP>_samples.txt
# Detecta si es continental (AFR/EUR/AMR/EAS/SAS) => usa columna super_pop (col 3)
# Si no => usa columna pop (col 2)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$PROJ_DIR/env.sh"

if [[ $# -ne 1 ]]; then
  echo "[X] Uso: $0 <POPS_FILE>" >&2
  exit 1
fi

POPS_FILE="$1"
[[ -s "$POPS_FILE" ]] || { echo "[X] No existe/está vacío: $POPS_FILE" >&2; exit 1; }
[[ -s "$PANEL_FILE" ]] || { echo "[X] No existe PANEL_FILE: $PANEL_FILE" >&2; exit 1; }

mapfile -t POPS < <(awk 'NF{gsub("\r",""); gsub(/^[ \t]+|[ \t]+$/,""); print toupper($0)}' "$POPS_FILE")

# Detecta si es continental
use_col=2
for p in "${POPS[@]}"; do
  if [[ "$p" =~ ^(AFR|EUR|AMR|EAS|SAS)$ ]]; then
    use_col=3
    break
  fi
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📋 Generando listas de muestras por población"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "[INFO] Panel: $PANEL_FILE"
echo "[INFO] Usando columna $use_col ($( [[ $use_col -eq 3 ]] && echo "super_pop" || echo "pop" ))"
echo

mkdir -p "$SAMP_DIR"

# Panel tiene header: sample pop super_pop gender ...
# sample=$1, pop=$2, super_pop=$3
for POP in "${POPS[@]}"; do
  out="$SAMP_DIR/${POP}_samples.txt"
  awk -v pop="$POP" -v col="$use_col" 'NR>1 && $col==pop {print $1}' "$PANEL_FILE" > "$out" || true
  n=$(wc -l < "$out" | tr -d ' ')
  if [[ "$n" -eq 0 ]]; then
    echo "  ❌ $POP: No encontrada"
    rm -f "$out"
  else
    echo "  ✓ $POP: $n muestras → $out"
  fi
done

echo
echo "✅ Listas listas en: $SAMP_DIR"

