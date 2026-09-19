#!/bin/bash
# ==============================================================================
# DESCARGA: 1000 Genomes Phase 3 - Genomas Completos (chr1-22)
# ==============================================================================
# Tamaño total: ~500GB
# Tiempo estimado: 1-2 días (depende de tu conexión)
# ==============================================================================

set -euo pipefail

# Configuración
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
DEST_DIR="$HOME/Documentos/Ancestria/data/vcf"
LOG_FILE="$HOME/Documentos/Ancestria/logs/download_1kgp.log"

mkdir -p "$DEST_DIR"
mkdir -p "$(dirname "$LOG_FILE")"

echo "=============================================="
echo "DESCARGA: 1000 Genomes Phase 3 (chr1-22)"
echo "Destino: $DEST_DIR"
echo "Log: $LOG_FILE"
echo "=============================================="
echo ""

# Función para descargar con reintentos
download_with_retry() {
    local url=$1
    local output=$2
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        echo "[$(date)] Intento $attempt de $max_attempts: $output" | tee -a "$LOG_FILE"
        
        if wget -c -O "$output" "$url" 2>&1 | tee -a "$LOG_FILE"; then
            echo "[$(date)] ✓ Descargado: $output" | tee -a "$LOG_FILE"
            return 0
        else
            echo "[$(date)] ✗ Fallo en intento $attempt" | tee -a "$LOG_FILE"
            attempt=$((attempt + 1))
            sleep 10
        fi
    done
    
    echo "[$(date)] ✗✗✗ ERROR: No se pudo descargar $output después de $max_attempts intentos" | tee -a "$LOG_FILE"
    return 1
}

# Verificar espacio disponible
available_space=$(df -BG "$DEST_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
required_space=200

if [ "$available_space" -lt "$required_space" ]; then
    echo "ERROR: Espacio insuficiente."
    echo "Disponible: ${available_space}GB"
    echo "Requerido: ${required_space}GB"
    exit 1
fi

echo "Espacio disponible: ${available_space}GB ✓"
echo ""

# Descargar chr1-22 + índices
for chr in {1..22}; do
    vcf_file="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    tbi_file="${vcf_file}.tbi"
    
    vcf_url="${BASE_URL}/${vcf_file}"
    tbi_url="${BASE_URL}/${tbi_file}"
    
    vcf_dest="${DEST_DIR}/${vcf_file}"
    tbi_dest="${DEST_DIR}/${tbi_file}"
    
    echo "=============================================="
    echo "CROMOSOMA $chr"
    echo "=============================================="
    
    # Descargar VCF
    if [ -f "$vcf_dest" ]; then
        echo "✓ Ya existe: $vcf_file"
    else
        download_with_retry "$vcf_url" "$vcf_dest" || exit 1
    fi
    
    # Descargar índice
    if [ -f "$tbi_dest" ]; then
        echo "✓ Ya existe: $tbi_file"
    else
        download_with_retry "$tbi_url" "$tbi_dest" || exit 1
    fi
    
    # Verificar integridad
    if [ -f "$vcf_dest" ] && [ -f "$tbi_dest" ]; then
        echo "✓ Cromosoma $chr completo"
        echo ""
    fi
done

echo ""
echo "=============================================="
echo "DESCARGA COMPLETADA"
echo "=============================================="
echo ""
echo "Archivos descargados:"
ls -lh "$DEST_DIR"/ALL.chr*.vcf.gz | wc -l
echo ""
echo "Espacio usado:"
du -sh "$DEST_DIR"
echo ""
echo "Siguiente paso: Construir paneles de referencia"
echo "  ./build_panels_10runs.sh"
