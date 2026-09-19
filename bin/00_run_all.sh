#!/usr/bin/env bash
# =============================================================================
# 00_run_all.sh
#
# Script maestro — corre el pipeline completo en orden.
# Puede retomarse desde cualquier punto si algo falla
# (todos los scripts verifican si la salida ya existe antes de recalcular).
#
# Tiempo estimado total: 2-3 días corriendo sin interrupciones.
#
# Uso: nohup bash 00_run_all.sh 2>&1 | tee logs/pipeline_$(date +%Y%m%d_%H%M).log &
# =============================================================================
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS=20
START_TIME=$(date +%s)

mkdir -p "$PROJ/logs"
LOG="$PROJ/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"

echo "============================================================" | tee -a "$LOG"
echo " Pipeline Ancestría Genética — Genoma completo"            | tee -a "$LOG"
echo " Inicio: $(date)"                                          | tee -a "$LOG"
echo " Threads: $THREADS"                                        | tee -a "$LOG"
echo " Log: $LOG"                                                | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

# Verificar dependencias
echo "" | tee -a "$LOG"
echo "[CHECK] Verificando dependencias..." | tee -a "$LOG"
for tool in gatk bcftools samtools CrossMap python3 wget; do
  if command -v "$tool" &>/dev/null; then
    echo "  ✓ $tool: $(command -v $tool)" | tee -a "$LOG"
  else
    echo "  ✗ $tool: NO ENCONTRADO" | tee -a "$LOG"
    exit 1
  fi
done

# Verificar Python packages
python3 -c "import numpy" 2>/dev/null && echo "  ✓ numpy" | tee -a "$LOG" || \
  { echo "  ✗ numpy — instala con: pip install numpy"; exit 1; }

# Verificar archivos de entrada
echo "" | tee -a "$LOG"
echo "[CHECK] Verificando archivos de entrada..." | tee -a "$LOG"
required_files=(
  "$PROJ/data/bam/HOSPITAL.bam"
  "$PROJ/data/bam/01-50_samples_UDB-102_482263.merged.bam"
  "$PROJ/data/bam/51-100_samples_UDB-103_482264.merged.bam"
  "$PROJ/data/ref/hg19.fa"
  "$PROJ/data/ref/hg38.fa"
  "$PROJ/data/ref/hg19ToHg38.over.chain.gz"
  "$PROJ/data/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz"
  "$PROJ/data/panel/gnomad_hgdp_1kg_meta.tsv"
)

all_ok=true
for f in "${required_files[@]}"; do
  if [[ -s "$f" ]]; then
    echo "  ✓ $(basename $f)" | tee -a "$LOG"
  else
    echo "  ✗ FALTA: $f" | tee -a "$LOG"
    all_ok=false
  fi
done

if [[ "$all_ok" == false ]]; then
  echo "" | tee -a "$LOG"
  echo "[!] Faltan archivos de entrada. Revisa los paths en env.sh" | tee -a "$LOG"
  echo "    Archivos de referencia necesarios:" | tee -a "$LOG"
  echo "    hg19.fa  : wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" | tee -a "$LOG"
  echo "    hg38.fa  : wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" | tee -a "$LOG"
  echo "    chain    : wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" | tee -a "$LOG"
  exit 1
fi

# Verificar espacio en disco
AVAIL_GB=$(df "$PROJ" | tail -1 | awk '{print int($4/1024/1024)}')
echo "" | tee -a "$LOG"
echo "[CHECK] Espacio disponible: ${AVAIL_GB} GB" | tee -a "$LOG"
if [[ "$AVAIL_GB" -lt 80 ]]; then
  echo "  [!] ADVERTENCIA: menos de 80 GB libres. El pipeline puede fallar." | tee -a "$LOG"
  echo "  Necesitas al menos 80 GB para procesar un VCF de gnomAD a la vez." | tee -a "$LOG"
fi

# ─────────────────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
echo " PASO 1/4 — GATK HaplotypeCaller (genoma completo)"   | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
T1=$(date +%s)

bash "$PROJ/bin/01_gatk_genome.sh" "$THREADS" 2>&1 | tee -a "$LOG"

T2=$(date +%s)
echo " Tiempo PASO 1: $(( (T2-T1)/3600 ))h $(( ((T2-T1)%3600)/60 ))min" | tee -a "$LOG"

# ─────────────────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
echo " PASO 2/4 — Panel gnomAD N5 (genoma completo)"         | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
T1=$(date +%s)

bash "$PROJ/bin/02_build_panel_genome.sh" "$THREADS" 2>&1 | tee -a "$LOG"

T2=$(date +%s)
echo " Tiempo PASO 2: $(( (T2-T1)/3600 ))h $(( ((T2-T1)%3600)/60 ))min" | tee -a "$LOG"

# ─────────────────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
echo " PASO 3/4 — iAdmix N5 (10 corridas)"                   | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
T1=$(date +%s)

bash "$PROJ/bin/03_iadmix_genome.sh" "$THREADS" 2>&1 | tee -a "$LOG"

T2=$(date +%s)
echo " Tiempo PASO 3: $(( (T2-T1)/3600 ))h $(( ((T2-T1)%3600)/60 ))min" | tee -a "$LOG"

# ─────────────────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
echo " PASO 4/4 — Métricas de distancia"                     | tee -a "$LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$LOG"
T1=$(date +%s)

bash "$PROJ/bin/04_distances_genome.sh" "$THREADS" 2>&1 | tee -a "$LOG"

T2=$(date +%s)
echo " Tiempo PASO 4: $(( (T2-T1)/3600 ))h $(( ((T2-T1)%3600)/60 ))min" | tee -a "$LOG"

# ─────────────────────────────────────────────────────────────────────────────
END_TIME=$(date +%s)
TOTAL=$(( END_TIME - START_TIME ))

echo "" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo " ✅ PIPELINE COMPLETADO"                                    | tee -a "$LOG"
echo " Fin: $(date)"                                             | tee -a "$LOG"
echo " Tiempo total: $(( TOTAL/3600 ))h $(( (TOTAL%3600)/60 ))min" | tee -a "$LOG"
echo ""                                                          | tee -a "$LOG"
echo " Resultados:"                                             | tee -a "$LOG"
echo "   iAdmix   : $PROJ/results/iadmix_genome/summary_N5_genome.tsv" | tee -a "$LOG"
echo "   Distancias: $PROJ/results/distances_genome/summary_distances_genome.tsv" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
