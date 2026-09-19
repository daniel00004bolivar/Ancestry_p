#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
VCF_DIR="$PROJ_DIR/work/gatk_genome"
PANEL_3WAY="$PROJ_DIR/refpanels/PANEL_3WAY_run1_extreme.txt"
PANEL_N11="$PROJ_DIR/refpanels/PANEL_N11_FULL_run1.txt"
OUTPUT_DIR="$PROJ_DIR/results/genetic_distances"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "CÁLCULO DE DISTANCIAS GENÉTICAS"
echo "=========================================="

# ============================================
# Extraer frecuencias de VCFs
# ============================================

echo "Extrayendo frecuencias de pools..."
for sample in POOL1 POOL2 HOSPITAL; do
    vcf="$VCF_DIR/${sample}_hg38.vcf.gz"
    out="$OUTPUT_DIR/${sample}_frequencies.txt"
    
    if [[ ! -s "$out" ]]; then
        echo "  $sample..."
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$vcf" | \
        awk 'BEGIN{OFS="\t"} {
            key = $1":"$2":"$3":"$4
            af = $5
            if(af >= 0 && af <= 1) print key, af
        }' > "$out"
    fi
done

# ============================================
# Panel 3-WAY
# ============================================

echo ""
echo "Preparando panel 3-WAY..."

# Verificar columnas del panel
head -1 "$PANEL_3WAY" | tr '\t' '\n' | nl

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1; i<=NF; i++) {
        if($i=="AFR") col_afr=i
        if($i=="EUR") col_eur=i  
        if($i=="AMR") col_amr=i
    }
    next
}
NR>1 && col_afr && col_eur && col_amr {
    key = $1":"$2":"$4":"$5
    print key, $col_afr, $col_eur, $col_amr
}' "$PANEL_3WAY" > "$OUTPUT_DIR/ref_3way_frequencies.txt"

# ============================================
# Panel N11
# ============================================

echo "Preparando panel N11..."

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    for(i=1; i<=NF; i++) {
        if($i=="IBS") c[1]=i
        if($i=="CEU") c[2]=i
        if($i=="YRI") c[3]=i
        if($i=="ESN") c[4]=i
        if($i=="GWD") c[5]=i
        if($i=="LWK") c[6]=i
        if($i=="PEL") c[7]=i
        if($i=="MXL") c[8]=i
        if($i=="PUR") c[9]=i
        if($i=="FIN") c[10]=i
        if($i=="CLM") c[11]=i
    }
    next
}
NR>1 {
    key = $1":"$2":"$4":"$5
    
    # Verificar que todas las columnas existan
    valid = 1
    for(i=1; i<=11; i++) {
        if(!c[i] || $c[i] == "") valid = 0
    }
    
    if(valid) {
        printf "%s", key
        for(i=1; i<=11; i++) printf "\t%s", $c[i]
        printf "\n"
    }
}' "$PANEL_N11" > "$OUTPUT_DIR/ref_N11_frequencies.txt"

echo "  Líneas procesadas: $(wc -l < $OUTPUT_DIR/ref_N11_frequencies.txt)"

# ============================================
# Script Python
# ============================================

cat > "$OUTPUT_DIR/calc_distances.py" << 'PYTHON'
import sys
import math

def calculate_metrics(freq_pool_file, freq_ref_file, pop_cols, output_file):
    pool_freqs = {}
    with open(freq_pool_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    pool_freqs[parts[0]] = float(parts[1])
                except:
                    continue
    
    ref_freqs = {}
    with open(freq_ref_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= len(pop_cols) + 1:
                key = parts[0]
                try:
                    freqs = [float(parts[i+1].strip()) for i in range(len(pop_cols))]
                    ref_freqs[key] = freqs
                except ValueError as e:
                    print(f"Error parsing line: {line.strip()}", file=sys.stderr)
                    print(f"  Parts: {parts[:5]}", file=sys.stderr)
                    continue
    
    common_snps = set(pool_freqs.keys()) & set(ref_freqs.keys())
    
    results = []
    for pop_idx, pop_name in enumerate(pop_cols):
        n = 0
        sum_ab = sum_ht = sum_hs = sum_nei_num = sum_bhatt = 0
        
        for snp in common_snps:
            pp = pool_freqs[snp]
            pr = ref_freqs[snp][pop_idx]
            
            if pr < 0 or pr > 1 or pp < 0 or pp > 1:
                continue
                
            n += 1
            sum_ab += pp * pr + (1 - pp) * (1 - pr)
            
            p_mean = (pp + pr) / 2
            sum_ht += 2 * p_mean * (1 - p_mean)
            sum_hs += pp * (1 - pp) + pr * (1 - pr)
            sum_nei_num += pp * pp + pr * pr
            sum_bhatt += math.sqrt(pp * pr) + math.sqrt((1 - pp) * (1 - pr))
        
        if n == 0:
            continue
            
        ab_mean = sum_ab / n
        fst = max(0, (sum_ht - sum_hs / 2) / sum_ht if sum_ht > 0 else 0)
        reynolds = -math.log(1 - fst) if fst < 1 else float('inf')
        nei_dist = -math.log(sum_nei_num / n) if sum_nei_num > 0 else 0
        bhatt_dist = sum_bhatt / n
        
        results.append({
            'pop': pop_name, 'ab': ab_mean, 'fst': fst,
            'reynolds': reynolds, 'nei': nei_dist,
            'bhatt': bhatt_dist, 'n': n
        })
    
    with open(output_file, 'w') as f:
        f.write("Pop\tÃ_B\tFST\tReynolds\tNei\tBhattacharyya\tn_SNPs\n")
        for r in results:
            f.write(f"{r['pop']}\t{r['ab']:.6f}\t{r['fst']:.6f}\t"
                   f"{r['reynolds']:.6f}\t{r['nei']:.6f}\t"
                   f"{r['bhatt']:.6f}\t{r['n']}\n")

if __name__ == "__main__":
    calculate_metrics(sys.argv[1], sys.argv[2], 
                     sys.argv[3].split(','), sys.argv[4])
PYTHON

# ============================================
# Calcular distancias
# ============================================

echo ""
echo "=== PANEL 3-WAY ===" | tee "$OUTPUT_DIR/summary_3way.txt"

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:" | tee -a "$OUTPUT_DIR/summary_3way.txt"
    python3 "$OUTPUT_DIR/calc_distances.py" \
        "$OUTPUT_DIR/${sample}_frequencies.txt" \
        "$OUTPUT_DIR/ref_3way_frequencies.txt" \
        "AFR,EUR,AMR" \
        "$OUTPUT_DIR/${sample}_distances_3way.txt"
    cat "$OUTPUT_DIR/${sample}_distances_3way.txt" | column -t | tee -a "$OUTPUT_DIR/summary_3way.txt"
    echo "" | tee -a "$OUTPUT_DIR/summary_3way.txt"
done

echo ""
echo "=== PANEL N11 ===" | tee "$OUTPUT_DIR/summary_N11.txt"

for sample in POOL1 POOL2 HOSPITAL; do
    echo "$sample:" | tee -a "$OUTPUT_DIR/summary_N11.txt"
    python3 "$OUTPUT_DIR/calc_distances.py" \
        "$OUTPUT_DIR/${sample}_frequencies.txt" \
        "$OUTPUT_DIR/ref_N11_frequencies.txt" \
        "IBS,CEU,YRI,ESN,GWD,LWK,PEL,MXL,PUR,FIN,CLM" \
        "$OUTPUT_DIR/${sample}_distances_N11.txt"
    cat "$OUTPUT_DIR/${sample}_distances_N11.txt" | column -t | tee -a "$OUTPUT_DIR/summary_N11.txt"
    echo "" | tee -a "$OUTPUT_DIR/summary_N11.txt"
done

echo "✓ COMPLETADO"
