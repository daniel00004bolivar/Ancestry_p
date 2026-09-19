#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
WORK="$PROJ_DIR/work"
VCF_DIR="$WORK/gatk_genome"
FREQ_DIR="$WORK/freq_genome"
RESULTS_DIR="$PROJ_DIR/results/genetic_distances"

mkdir -p "$RESULTS_DIR" "$FREQ_DIR"

declare -A SAMPLES=(
  [HOSPITAL]="$VCF_DIR/HOSPITAL_hg38.vcf.gz"
  [POOL1]="$VCF_DIR/POOL1_hg38.vcf.gz"
  [POOL2]="$VCF_DIR/POOL2_hg38.vcf.gz"
)

echo "============================================================"
echo " Distancias genéticas — Genoma completo"
echo " N3: AFR, EUR, AMR"
echo " N5: AFR, EUR, AMR, EAS, SAS"
echo " N11: IBS, CEU, YRI, ESN, GWD, LWK, PEL, MXL, PUR, FIN, CLM"
echo "============================================================"

# ============================================
# Script Python para cálculos
# ============================================

cat > "$RESULTS_DIR/calc_distances.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import numpy as np
import sys

def load_freq(path):
    """Carga frecuencias en diccionario {chr:pos:ref:alt -> af}"""
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
    """Calcula Ã_B, FST, Reynolds, Nei, Bhattacharyya"""
    n = len(p_pool)
    if n == 0:
        return None

    eps = 1e-9

    # Ã_B (Weir & Goudet 2026)
    ab = np.mean(p_pool * p_ref + (1 - p_pool) * (1 - p_ref))

    # FST de Weir-Cockerham
    p_bar = (p_pool + p_ref) / 2.0
    H_T = 2 * p_bar * (1 - p_bar)
    H_S = (2*p_pool*(1-p_pool) + 2*p_ref*(1-p_ref)) / 2.0
    
    valid = H_T > eps
    fst = np.sum((H_T - H_S)[valid]) / np.sum(H_T[valid])
    fst = max(0.0, fst)

    # Reynolds
    reynolds = -np.log(max(1 - fst, eps))

    # Nei (1972)
    J12 = np.mean(p_pool*p_ref + (1-p_pool)*(1-p_ref))
    J1  = np.mean(p_pool**2 + (1-p_pool)**2)
    J2  = np.mean(p_ref**2 + (1-p_ref)**2)
    denom_nei = np.sqrt(J1 * J2)
    nei = -np.log(J12 / denom_nei) if denom_nei > eps else np.nan

    # Bhattacharyya
    bc = np.mean(np.sqrt(p_pool * p_ref) + np.sqrt((1-p_pool) * (1-p_ref)))
    bhatt = -np.log(max(bc / 2.0, eps))

    return {
        'n_snps': n,
        'AB': round(ab, 6),
        'FST': round(fst, 6),
        'Reynolds': round(reynolds, 6),
        'Nei': round(nei, 6),
        'Bhattacharyya': round(bhatt, 6),
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pool_freq', required=True)
    parser.add_argument('--ref_freq', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--pop', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    pool_d = load_freq(args.pool_freq)
    ref_d = load_freq(args.ref_freq)

    common_keys = set(pool_d.keys()) & set(ref_d.keys())
    
    if len(common_keys) == 0:
        sys.stderr.write(f"[!] Sin SNPs comunes: {args.sample} vs {args.pop}\n")
        return

    p_pool = np.array([pool_d[k] for k in common_keys], dtype=np.float64)
    p_ref = np.array([ref_d[k] for k in common_keys], dtype=np.float64)

    p_pool = np.clip(p_pool, 1e-6, 1-1e-6)
    p_ref = np.clip(p_ref, 1e-6, 1-1e-6)

    metrics = calc_metrics(p_pool, p_ref)
    if metrics is None:
        return

    with open(args.out, 'w') as f:
        f.write("Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya\n")
        f.write(f"{args.sample}\t{args.pop}\t{metrics['n_snps']}\t"
                f"{metrics['AB']}\t{metrics['FST']}\t{metrics['Reynolds']}\t"
                f"{metrics['Nei']}\t{metrics['Bhattacharyya']}\n")

    sys.stderr.write(f"[✓] {args.sample} vs {args.pop}: FST={metrics['FST']:.4f} AB={metrics['AB']:.4f} (n={metrics['n_snps']})\n")

if __name__ == '__main__':
    main()
PYEOF

chmod +x "$RESULTS_DIR/calc_distances.py"

# ============================================
# [1] Extraer frecuencias de pools
# ============================================

echo ""
echo "[1/4] Extrayendo frecuencias de pools..."

for sample in "${!SAMPLES[@]}"; do
    vcf="${SAMPLES[$sample]}"
    out_tsv="$RESULTS_DIR/${sample}_pool_af.tsv"

    if [[ -s "$out_tsv" ]]; then
        echo "  [=] $sample ya extraído"
        continue
    fi

    echo "  [...] $sample"
    bcftools query \
        -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
        -i 'TYPE="snp" && INFO/AF>0' \
        "$vcf" | \
    awk 'BEGIN{OFS="\t"} NF==5 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/ && $3!=$4 {print}' > "$out_tsv"
    
    echo "  [✓] $sample: $(wc -l < $out_tsv) SNPs"
done

# ============================================
# [2] Extraer frecuencias de referencia
# ============================================

echo ""
echo "[2/4] Extrayendo frecuencias de paneles de referencia..."

# Panel N3
for pop in AFR EUR AMR; do
    out="$FREQ_DIR/${pop}_n3_af.tsv"
    if [[ ! -s "$out" ]]; then
        echo "  [...] $pop (N3)"
        awk -v pop="$pop" 'BEGIN{FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) if($i==pop) col=i; next}
        col && $col!="" {print $1,$2,$4,$5,$col}' \
        "$PROJ_DIR/refpanels/PANEL_3WAY_run1_extreme.txt" > "$out"
        echo "  [✓] $pop: $(wc -l < $out) SNPs"
    fi
done

# Panel N5
for pop in AFR EUR AMR EAS SAS; do
    out="$FREQ_DIR/${pop}_n5_af.tsv"
    if [[ ! -s "$out" ]]; then
        echo "  [...] $pop (N5)"
        awk -v pop="$pop" 'BEGIN{FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) if($i==pop) col=i; next}
        col && $col!="" {print $1,$2,$4,$5,$col}' \
        "$PROJ_DIR/refpanels/FINAL_FREQUENCIESN5.txt" > "$out"
        echo "  [✓] $pop: $(wc -l < $out) SNPs"
    fi
done

# Panel N11
for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
    out="$FREQ_DIR/${pop}_n11_af.tsv"
    if [[ ! -s "$out" ]]; then
        echo "  [...] $pop (N11)"
        awk -v pop="$pop" 'BEGIN{FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) if($i==pop) col=i; next}
        col && $col!="" {print $1,$2,$4,$5,$col}' \
        "$PROJ_DIR/refpanels/FINAL_FREQUENCIESN11.txt" > "$out"
        echo "  [✓] $pop: $(wc -l < $out) SNPs"
    fi
done

# ============================================
# [3] Calcular distancias (paralelo)
# ============================================

echo ""
echo "[3/4] Calculando distancias..."

pids=()

# N3
echo "  Panel N3..."
for sample in "${!SAMPLES[@]}"; do
    for pop in AFR EUR AMR; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n3.tsv"
        [[ -s "$out" ]] && continue
        
        python3 "$RESULTS_DIR/calc_distances.py" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$FREQ_DIR/${pop}_n3_af.tsv" \
            --sample "$sample" --pop "$pop" --out "$out" &
        pids+=($!)
        
        if [[ ${#pids[@]} -ge 10 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

# N5
echo "  Panel N5..."
for sample in "${!SAMPLES[@]}"; do
    for pop in AFR EUR AMR EAS SAS; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n5.tsv"
        [[ -s "$out" ]] && continue
        
        python3 "$RESULTS_DIR/calc_distances.py" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$FREQ_DIR/${pop}_n5_af.tsv" \
            --sample "$sample" --pop "$pop" --out "$out" &
        pids+=($!)
        
        if [[ ${#pids[@]} -ge 10 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

# N11
echo "  Panel N11..."
for sample in "${!SAMPLES[@]}"; do
    for pop in IBS CEU YRI ESN GWD LWK PEL MXL PUR FIN CLM; do
        out="$RESULTS_DIR/${sample}_vs_${pop}_n11.tsv"
        [[ -s "$out" ]] && continue
        
        python3 "$RESULTS_DIR/calc_distances.py" \
            --pool_freq "$RESULTS_DIR/${sample}_pool_af.tsv" \
            --ref_freq "$FREQ_DIR/${pop}_n11_af.tsv" \
            --sample "$sample" --pop "$pop" --out "$out" &
        pids+=($!)
        
        if [[ ${#pids[@]} -ge 10 ]]; then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
done

for pid in "${pids[@]}"; do wait "$pid"; done

# ============================================
# [4] Consolidar resultados
# ============================================

echo ""
echo "[4/4] Consolidando resultados..."

for panel in n3 n5 n11; do
    summary="$RESULTS_DIR/summary_${panel}.tsv"
    echo -e "Sample\tPop\tn_SNPs\tAB\tFST\tReynolds\tNei\tBhattacharyya" > "$summary"
    
    for f in "$RESULTS_DIR"/*_vs_*_${panel}.tsv; do
        [[ -s "$f" ]] && tail -1 "$f" >> "$summary"
    done
    
    echo ""
    echo "=== PANEL ${panel^^} ===" 
    cat "$summary" | column -t
done

echo ""
echo "============================================================"
echo " ✅ Distancias completadas"
echo "============================================================"
ls -lh "$RESULTS_DIR"/summary_*.tsv
