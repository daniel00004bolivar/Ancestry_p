#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/work/gatk_genome"
OUTPUT_DIR="$PROJ_DIR/results/genetic_distances"

# Paneles
PANEL_3WAY="$PROJ_DIR/refpanels/PANEL_3WAY_run1_extreme.txt"
PANEL_N5="$PROJ_DIR/refpanels/FINAL_FREQUENCIESN5.txt"
PANEL_N11="$PROJ_DIR/refpanels/FINAL_FREQUENCIESN11.txt"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "DISTANCIAS GENÉTICAS - N3, N5, N11"
echo "=========================================="

# ============================================
# Script Python para cálculos
# ============================================

cat > "$OUTPUT_DIR/calc.py" << 'PYTHON'
import sys, math

def calc(pool_file, ref_file, pops, out_file):
    pool = {}
    with open(pool_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 2:
                try: pool[p[0]] = float(p[1])
                except: pass
    
    ref = {}
    with open(ref_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= len(pops) + 1:
                try:
                    ref[p[0]] = [float(p[i+1]) for i in range(len(pops))]
                except: pass
    
    common = set(pool.keys()) & set(ref.keys())
    results = []
    
    for idx, pop in enumerate(pops):
        n = sum_ab = sum_ht = sum_hs = sum_nei = sum_bhatt = 0
        
        for snp in common:
            pp, pr = pool[snp], ref[snp][idx]
            if not (0 <= pp <= 1 and 0 <= pr <= 1): continue
            
            n += 1
            sum_ab += pp*pr + (1-pp)*(1-pr)
            pm = (pp+pr)/2
            sum_ht += 2*pm*(1-pm)
            sum_hs += pp*(1-pp) + pr*(1-pr)
            sum_nei += pp*pp + pr*pr
            sum_bhatt += math.sqrt(pp*pr) + math.sqrt((1-pp)*(1-pr))
        
        if n == 0: continue
        
        fst = max(0, (sum_ht - sum_hs/2) / sum_ht if sum_ht > 0 else 0)
        results.append({
            'pop': pop,
            'ab': sum_ab/n,
            'fst': fst,
            'rey': -math.log(1-fst) if fst < 1 else 999,
            'nei': -math.log(sum_nei/n),
            'bhat': sum_bhatt/n,
            'n': n
        })
    
    with open(out_file, 'w') as f:
        f.write("Pop\tÃ_B\tFST\tReynolds\tNei\tBhatt\tn\n")
        for r in results:
            f.write(f"{r['pop']}\t{r['ab']:.6f}\t{r['fst']:.6f}\t{r['rey']:.6f}\t{r['nei']:.6f}\t{r['bhat']:.6f}\t{r['n']}\n")

calc(sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4])
PYTHON

# ============================================
# Extraer frecuencias de pools
# ============================================

echo "Extrayendo frecuencias..."
for sample in POOL1 POOL2 HOSPITAL; do
    out="$OUTPUT_DIR/${sample}_freq.txt"
    [[ -s "$out" ]] && continue
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$VCF_DIR/${sample}_hg38.vcf.gz" | \
    awk 'BEGIN{OFS="\t"} $5>=0 && $5<=1 {print $1":"$2":"$3":"$4, $5}' > "$out"
done

# ============================================
# Preparar paneles de referencia
# ============================================

# N3
awk 'NR==1{for(i=1;i<=NF;i++){if($i=="AFR")a=i;if($i=="EUR")e=i;if($i=="AMR")m=i}next}
     {print $1":"$2":"$4":"$5"\t"$a"\t"$e"\t"$m}' "$PANEL_3WAY" > "$OUTPUT_DIR/ref_n3.txt"

# N5  
awk 'NR==1{for(i=1;i<=NF;i++){if($i=="AFR")a=i;if($i=="NFE")n=i;if($i=="AMR")m=i;if($i=="EAS")ea=i;if($i=="SAS")s=i}next}
     {print $1":"$2":"$4":"$5"\t"$a"\t"$n"\t"$m"\t"$ea"\t"$s}' "$PANEL_N5" > "$OUTPUT_DIR/ref_n5.txt"

# N11
awk 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
     {printf "%s:%s:%s:%s", $1,$2,$4,$5;
      for(p in c) if(p~/IBS|CEU|YRI|ESN|GWD|LWK|PEL|MXL|PUR|FIN|CLM/) printf "\t%s",$c[p];
      print ""}' "$PANEL_N11" > "$OUTPUT_DIR/ref_n11.txt"

echo "Paneles preparados"

# ============================================
# Calcular distancias
# ============================================

for panel in n3 n5 n11; do
    [[ $panel == "n3" ]] && pops="AFR,EUR,AMR"
    [[ $panel == "n5" ]] && pops="AFR,NFE,AMR,EAS,SAS"
    [[ $panel == "n11" ]] && pops="IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM"
    
    echo ""
    echo "=== PANEL ${panel^^} ===" | tee "$OUTPUT_DIR/summary_${panel}.txt"
    
    for sample in POOL1 POOL2 HOSPITAL; do
        echo "$sample:" | tee -a "$OUTPUT_DIR/summary_${panel}.txt"
        python3 "$OUTPUT_DIR/calc.py" \
            "$OUTPUT_DIR/${sample}_freq.txt" \
            "$OUTPUT_DIR/ref_${panel}.txt" \
            "$pops" \
            "$OUTPUT_DIR/${sample}_${panel}.txt"
        cat "$OUTPUT_DIR/${sample}_${panel}.txt" | column -t | tee -a "$OUTPUT_DIR/summary_${panel}.txt"
        echo "" | tee -a "$OUTPUT_DIR/summary_${panel}.txt"
    done
done

echo "✓ COMPLETADO"
ls -lh "$OUTPUT_DIR"/summary_*.txt
