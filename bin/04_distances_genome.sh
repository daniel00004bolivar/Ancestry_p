#!/usr/bin/env bash
# =============================================================================
# 04_distances_genome.sh
#
# Calcula Ã_B, FST, Reynolds, Nei y Bhattacharyya para genoma completo.
# Usa los VCFs filtrados y liftoveados de 01_gatk_genome.sh
# y las frecuencias de referencia calculadas en 02_build_panel_genome.sh
#
# Uso: bash 04_distances_genome.sh [THREADS]
# =============================================================================
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$PROJ/env.sh"

THREADS="${1:-20}"
POPS=(AFR NFE AMR EAS SAS)
SAMP_DIR="$PROJ/data/samples/GNOMAD_GENOME"
FREQ_DIR="$WORK/freq_genome"
GATK_DIR="$WORK/gatk_genome"
RESULTS_DIR="$PROJ/results/distances_genome"
mkdir -p "$RESULTS_DIR"

declare -A SAMPLES=(
  [HOSPITAL]="$GATK_DIR/HOSPITAL_hg38.vcf.gz"
  [POOL1]="$GATK_DIR/POOL1_hg38.vcf.gz"
  [POOL2]="$GATK_DIR/POOL2_hg38.vcf.gz"
)

echo "============================================================"
echo " Métricas de distancia — Genoma completo"
echo " Muestras : ${!SAMPLES[*]}"
echo " Pobs ref : ${POPS[*]}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Para cada muestra y cada población de referencia:
# extraer frecuencias del pool y de la referencia, calcular métricas
# -----------------------------------------------------------------------------

# Script Python embebido para calcular todas las métricas
# Más eficiente que bash para operaciones vectorizadas sobre millones de SNPs
CALC_SCRIPT="$WORK/calc_distances.py"

cat > "$CALC_SCRIPT" << 'PYEOF'
#!/usr/bin/env python3
"""
Calcula Ã_B, FST (Weir-Cockerham), Reynolds, Nei, Bhattacharyya
entre un pool (frecuencias del VCF de GATK) y una población de referencia
(frecuencias calculadas desde genotipos individuales).

Entrada:
  --pool_freq   : TSV con columnas chr pos ref alt af_pool
  --ref_freq    : TSV con columnas chr pos ref alt af_ref
  --sample      : nombre de la muestra
  --pop         : nombre de la población de referencia
  --out         : archivo de salida con resultados

Salida:
  TSV con: sample, pop, n_snps, AB, FST, Reynolds, Nei, Bhattacharyya
"""
import argparse
import numpy as np
import sys

def load_freq(path):
    """Carga frecuencias en un diccionario {chr:pos:ref:alt -> af}"""
    freq = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            key = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"
            try:
                freq[key] = float(parts[4])
            except ValueError:
                continue
    return freq

def calc_metrics(p_pool, p_ref):
    """
    Calcula todas las métricas sobre arrays numpy de frecuencias.
    p_pool, p_ref : arrays 1D de frecuencias alélicas del alelo alternativo.
    """
    n = len(p_pool)
    if n == 0:
        return None

    eps = 1e-9  # evitar log(0)

    # ── Ã_B (Weir & Goudet 2026, Box 2) ──────────────────────────────────
    # Probabilidad de identidad alélica en estado entre pool y referencia
    ab = np.mean(p_pool * p_ref + (1 - p_pool) * (1 - p_ref))

    # ── FST de Weir-Cockerham ─────────────────────────────────────────────
    # p_bar = frecuencia media entre las dos "poblaciones"
    p_bar = (p_pool + p_ref) / 2.0

    # Heterocigosidades
    H_T = 2 * p_bar * (1 - p_bar)          # heterocigosidad total
    H_S = (2*p_pool*(1-p_pool) +            # promedio dentro de grupos
           2*p_ref *(1-p_ref )) / 2.0

    # FST global: suma de numeradores / suma de denominadores (no promedio de FSTs)
    num = H_T - H_S
    denom = H_T
    valid = denom > eps
    fst = np.sum(num[valid]) / np.sum(denom[valid])
    fst = max(0.0, fst)  # FST no puede ser negativo

    # ── Reynolds ─────────────────────────────────────────────────────────
    reynolds = -np.log(max(1 - fst, eps))

    # ── Nei (1972) ───────────────────────────────────────────────────────
    # J12 = probabilidad de identidad entre una pop y la otra
    J12 = np.mean(p_pool*p_ref + (1-p_pool)*(1-p_ref))
    J1  = np.mean(p_pool**2    + (1-p_pool)**2)
    J2  = np.mean(p_ref**2     + (1-p_ref)**2)
    denom_nei = np.sqrt(J1 * J2)
    nei = -np.log(J12 / denom_nei) if denom_nei > eps else np.nan

    # ── Bhattacharyya ────────────────────────────────────────────────────
    # BC = coeficiente de Bhattacharyya sobre distribución binaria
    bc = np.mean(
        np.sqrt(p_pool * p_ref) +
        np.sqrt((1-p_pool) * (1-p_ref))
    )
    bhatt = -np.log(max(bc / 2.0, eps))

    return {
        'n_snps': n,
        'AB':     round(ab,      6),
        'FST':    round(fst,     6),
        'Reynolds': round(reynolds, 6),
        'Nei':    round(nei,     6),
        'Bhattacharyya': round(bhatt, 6),
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pool_freq', required=True)
    parser.add_argument('--ref_freq',  required=True)
    parser.add_argument('--sample',    required=True)
    parser.add_argument('--pop',       required=True)
    parser.add_argument('--out',       required=True)
    args = parser.parse_args()

    sys.stderr.write(f"  Cargando {args.sample} vs {args.pop}...\n")

    pool_d = load_freq(args.pool_freq)
    ref_d  = load_freq(args.ref_freq)

    # Intersección de SNPs presentes en ambas fuentes
    common_keys = set(pool_d.keys()) & set(ref_d.keys())
    sys.stderr.write(f"  SNPs en común: {len(common_keys)}\n")

    if len(common_keys) == 0:
        sys.stderr.write(f"  [!] Sin SNPs en común para {args.sample} vs {args.pop}\n")
        return

    p_pool = np.array([pool_d[k] for k in common_keys], dtype=np.float64)
    p_ref  = np.array([ref_d[k]  for k in common_keys], dtype=np.float64)

    # Clamping por si acaso
    p_pool = np.clip(p_pool, 1e-6, 1-1e-6)
    p_ref  = np.clip(p_ref,  1e-6, 1-1e-6)

    metrics = calc_metrics(p_pool, p_ref)
    if metrics is None:
        return

    with open(args.out, 'w') as f:
        f.write("Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya\n")
        f.write(f"{args.sample}\t{args.pop}\t{metrics['n_snps']}\t"
                f"{metrics['AB']}\t{metrics['FST']}\t{metrics['Reynolds']}\t"
                f"{metrics['Nei']}\t{metrics['Bhattacharyya']}\n")

    sys.stderr.write(f"  [✓] {args.sample} vs {args.pop} — "
                     f"FST={metrics['FST']:.4f} AB={metrics['AB']:.4f}\n")

if __name__ == '__main__':
    main()
PYEOF

chmod +x "$CALC_SCRIPT"

# -----------------------------------------------------------------------------
# Extraer frecuencias del pool desde el VCF de GATK (campo AF del INFO)
# -----------------------------------------------------------------------------
echo ""
echo "[1/3] Extrayendo frecuencias del pool desde VCFs de GATK..."

for sample in "${!SAMPLES[@]}"; do
  vcf="${SAMPLES[$sample]}"
  out_tsv="$RESULTS_DIR/${sample}_pool_af.tsv"

  if [[ -s "$out_tsv" ]]; then
    echo "  [=] $sample pool AF ya extraído."
    continue
  fi

  if [[ ! -s "$vcf" ]]; then
    echo "  [!] VCF no encontrado para $sample: $vcf" >&2
    continue
  fi

  echo "  [...] Extrayendo AF del pool para $sample..."

  bcftools query \
    --threads "$THREADS" \
    -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
    -i 'TYPE="snp" && INFO/AF>0' \
    "$vcf" \
    | awk 'BEGIN{OFS="\t"}
        NF==5 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ && $3!=$4 {
          print $1,$2,$3,$4,$5
        }' > "$out_tsv"

  echo "  [✓] $sample: $(wc -l < "$out_tsv") SNPs"
done

# -----------------------------------------------------------------------------
# Combinar frecuencias de referencia por cromosoma en un TSV único por pop
# (los TSVs individuales ya existen de 02_build_panel_genome.sh)
# -----------------------------------------------------------------------------
echo ""
echo "[2/3] Consolidando frecuencias de referencia..."

for pop in "${POPS[@]}"; do
  out_tsv="$RESULTS_DIR/${pop}_ref_af_genome.tsv"

  if [[ -s "$out_tsv" ]]; then
    echo "  [=] $pop ref AF ya consolidado."
    continue
  fi

  echo "  [...] Consolidando $pop..."
  cat "$FREQ_DIR/${pop}_chr"*"_af.tsv" | sort -k1,1 -k2,2n > "$out_tsv"
  echo "  [✓] $pop: $(wc -l < "$out_tsv") SNPs"
done

# -----------------------------------------------------------------------------
# Calcular métricas: cada muestra × cada población (paralelo)
# -----------------------------------------------------------------------------
echo ""
echo "[3/3] Calculando métricas de distancia..."

pids=()
for sample in "${!SAMPLES[@]}"; do
  pool_tsv="$RESULTS_DIR/${sample}_pool_af.tsv"
  [[ ! -s "$pool_tsv" ]] && continue

  for pop in "${POPS[@]}"; do
    ref_tsv="$RESULTS_DIR/${pop}_ref_af_genome.tsv"
    out_file="$RESULTS_DIR/${sample}_vs_${pop}_metrics.tsv"

    if [[ -s "$out_file" ]]; then
      echo "  [=] $sample vs $pop ya calculado."
      continue
    fi

    python3 "$CALC_SCRIPT" \
      --pool_freq "$pool_tsv" \
      --ref_freq  "$ref_tsv" \
      --sample    "$sample" \
      --pop       "$pop" \
      --out       "$out_file" &

    pids+=($!)

    # Controlar paralelismo — máximo 15 procesos simultáneos
    if [[ ${#pids[@]} -ge 15 ]]; then
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
    fi
  done
done

for pid in "${pids[@]}"; do wait "$pid"; done

# -----------------------------------------------------------------------------
# Consolidar todos los resultados en una tabla resumen
# -----------------------------------------------------------------------------
SUMMARY="$RESULTS_DIR/summary_distances_genome.tsv"
echo "Sample	Pop	n_SNPs	AB	FST	Reynolds	Nei	Bhattacharyya" > "$SUMMARY"

for sample in "${!SAMPLES[@]}"; do
  for pop in "${POPS[@]}"; do
    f="$RESULTS_DIR/${sample}_vs_${pop}_metrics.tsv"
    [[ -s "$f" ]] && tail -1 "$f" >> "$SUMMARY"
  done
done

echo ""
echo "============================================================"
echo " ✅ Distancias completadas"
echo " Resumen: $SUMMARY"
echo "============================================================"
echo ""
cat "$SUMMARY" | column -t
