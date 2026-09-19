#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/work/gatk_genome"
FREQ_DIR="$PROJ_DIR/work/freq_genome"
RESULTS_DIR="$PROJ_DIR/results/genetic_distances"

mkdir -p "$RESULTS_DIR" "$FREQ_DIR"

declare -A SAMPLES=(
  [HOSPITAL]="$VCF_DIR/HOSPITAL_hg38.vcf.gz"
  [POOL1]="$VCF_DIR/POOL1_hg38.vcf.gz"
  [POOL2]="$VCF_DIR/POOL2_hg38.vcf.gz"
)

echo "============================================================"
echo " Distancias genéticas"
echo " N3: AFR, EUR, AMR"
echo " N11: IBS, CEU, YRI, ESN, GWD, LWK, PEL, MXL, PUR, FIN, CLM"
echo "============================================================"

# ============================================
# Script Python
# ============================================

cat > "$RESULTS_DIR/calc.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse, numpy as np, sys

def load_freq(path):
    freq = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 5: continue
            key = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"
            try: freq[key] = float(parts[4])
            except: continue
    return freq

def calc_metrics(p_pool, p_ref):
    n = len(p_pool)
    if n == 0: return None
    
    eps = 1e-9
    ab = np.mean(p_pool * p_ref + (1 - p_pool) * (1 - p_ref))
    
    p_bar = (p_pool + p_ref) / 2.0
    H_T = 2 * p_bar * (1 - p_bar)
    H_S = (2*p_pool*(1-p_pool) + 2*p_ref*(1-p_ref)) / 2.0
    valid = H_T > eps
    fst = np.sum((H_T - H_S)[valid]) / np.sum(H_T[valid])
    fst = max(0.0, fst)
    
    reynolds = -np.log(max(1 - fst, eps))
    
    J12 = np.mean(p_pool*p_ref + (1-p_pool)*(1-p_ref))
    J1 = np.mean(p_pool**2 + (1-p_pool)**2)
    J2 = np.mean(p_ref**2 + (1-p_ref)**2)
    denom_nei = np.sqrt(J1 * J2)
    nei = -np.log(J12 / denom_nei) if denom_nei > eps else np.nan
    
    bc = np.mean(np.sqrt(p_pool * p_ref) + np.sqrt((1-p_pool) * (1-p_ref)))
    bhatt = -np.log(max(bc / 2.0, eps))
    
    return {'n': n, 'AB': round(ab,6), 'FST': round(fst,6), 
            'Reynolds': round(reynolds,6), 'Nei': round(nei,6), 'Bhatt': round(bhatt,6)}

parser = argparse.ArgumentParser()
parser.add_argument('--pool', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--pop', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

pool_d = load_freq(args.pool)
ref_d = load_freq(args.ref)
common = set(pool_d.keys()) & set(ref_d.keys())

if len(common) == 0:
    sys.stderr.write(f"[!] Sin SNPs: {args.sample} vs {args.pop}\n")
    sys.exit(0)

p_pool = np.clip(np.array([pool_d[k] for k in common], dtype=np.float64), 1e-6, 1-1e-6)
p_ref = np.clip(np.array([ref_d[k] for k in common], dtype=np.float64), 1e-6, 1-1e-6)

m = calc_metrics(p_pool, p_ref)
if m:
    with open(args.out, 'w') as f:
        f.write("Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya\n")
        f.write(f"{args.sample}\t{args.pop}\t{m['n']}\t{m['AB']}\t{m['FST']}\t{m['Reynolds']}\t{m['Nei']}\t{m['Bhatt']}\n")
    sys.stderr.write(f"[✓] {args.sample} vs {args.pop}: FST={m['FST']:.4f} n={m['n']}\n")
PYEOF

chmod +x "$RESULTS_DIR/calc.py"

# ============================================
# [1] Extraer frecuencias pools
# ============================================

echo ""
echo "[1/3] Extrayendo frecuencias pools..."

for sample in "${!SAMPLES[@]}"; do
    out="$RESULTS_DIR/${sample}_pool.tsv"
    [[ -s "$out" ]] && continue
    
    echo "  $sample..."
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" -i 'TYPE="snp" && INFO/AF>0' "${SAMPLES[$sample]}" | \
    awk 'BEGIN{OFS="\t"} NF==5 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ && $3!=$4' > "$out"
done

# ============================================
# [2] Extraer frecuencias referencia
# ============================================

echo ""
echo "[2/3] Extrayendo frecuencias referencia..."

# N3
for pop in AFR EUR AMR; do
    out="$FREQ_DIR/${pop}_n3.tsv"
    [[ -s "$out" ]] && continue
    awk -v p="$pop" 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)if($i==p)c=i;next} c{print $1,$2,$4,$5,$c}' \
        "$PROJ_DIR/refpanels/PANEL_3WAY_run1_extreme.txt" > "$out"
done

# N11
for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
    out="$FREQ_DIR/${pop}_n11.tsv"
    [[ -s "$out" ]] && continue
    awk -v p="$pop" 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++)if($i==p)c=i;next} c{print $1,$2,$4,$5,$c}' \
        "$PROJ_DIR/refpanels/FINAL_FREQUENCIESN11.txt" > "$out"
done

# ============================================
# [3] Calcular distancias
# ============================================

echo ""
echo "[3/3] Calculando distancias..."

pids=()

# N3
for sample in "${!SAMPLES[@]}"; do
    for pop in AFR EUR AMR; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n3.tsv"
        [[ -s "$out" ]] && continue
        python3 "$RESULTS_DIR/calc.py" --pool "$RESULTS_DIR/${sample}_pool.tsv" \
            --ref "$FREQ_DIR/${pop}_n3.tsv" --sample "$sample" --pop "$pop" --out "$out" &
        pids+=($!)
        [[ ${#pids[@]} -ge 8 ]] && wait "${pids[0]}" && pids=("${pids[@]:1}")
    done
done

# N11
for sample in "${!SAMPLES[@]}"; do
    for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n11.tsv"
        [[ -s "$out" ]] && continue
        python3 "$RESULTS_DIR/calc.py" --pool "$RESULTS_DIR/${sample}_pool.tsv" \
            --ref "$FREQ_DIR/${pop}_n11.tsv" --sample "$sample" --pop "$pop" --out "$out" &
        pids+=($!)
        [[ ${#pids[@]} -ge 8 ]] && wait "${pids[0]}" && pids=("${pids[@]:1}")
    done
done

for pid in "${pids[@]}"; do wait "$pid"; done

# ============================================
# [4] Consolidar
# ============================================

echo ""
echo "Consolidando..."

for panel in n3 n11; do
    summ="$RESULTS_DIR/summary_${panel}.tsv"
    echo -e "Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya" > "$summ"
    for f in "$RESULTS_DIR"/*_vs_*_${panel}.tsv; do
        [[ -s "$f" ]] && tail -1 "$f" >> "$summ"
    done
    echo ""
    echo "=== ${panel^^} ==="
    cat "$summ" | column -t
done

echo ""
echo "✅ Completado"
