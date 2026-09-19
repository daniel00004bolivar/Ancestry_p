# PIPELINE COMPLETO: PANELES DE ANCESTRÍA CON 10 CORRIDAS

## 📋 RESUMEN

Este pipeline construye paneles de referencia científicamente rigurosos para análisis de ancestría:

- **Panel A: 3 Superpoblaciones** (AFR, EUR, AMR) × 10 corridas
- **Panel B: N11 (11 poblaciones)** (IBS, CEU, YRI, ESN, GWD, LWK, PEL, MXL, PUR, FIN, CLM) × 10 corridas

Cada corrida selecciona **100 individuos al azar CON reemplazo** de cada población, garantizando variabilidad estadística.

---

## ⏱️ TIEMPO TOTAL ESTIMADO

- **Descarga 1KGP Phase 3**: 1-2 días (500GB)
- **Construcción Panel A**: 6-8 horas
- **Construcción Panel B**: 8-10 horas
- **Ejecución iAdmix** (60 análisis): 1-2 días

**TOTAL: ~4-5 días**

---

## 📦 PREREQUISITOS

### Software necesario:
```bash
# Verificar instalaciones
bcftools --version  # >= 1.10
python2 --version   # Para iAdmix
wget --version
```

### Espacio en disco:
- VCFs: 500GB
- Paneles: 20GB
- Trabajo temporal: 100GB
- **TOTAL: ~650GB libres**

---

## 🚀 PASO 1: DESCARGAR 1000 GENOMES PHASE 3

```bash
cd ~/Documentos/Ancestria

# Copiar scripts al proyecto
cp download_1kgp_phase3.sh .
cp build_panels_master.sh .
cp build_panel_3way_10runs.sh .
cp build_panel_N11_10runs.sh .

# Dar permisos
chmod +x *.sh

# INICIAR DESCARGA (1-2 días)
./download_1kgp_phase3.sh
```

**Características:**
- Descarga chr1-22 + índices
- Reintentos automáticos si falla
- Log en `logs/download_1kgp.log`
- Verifica espacio disponible antes de empezar

**Monitorear progreso:**
```bash
# Ver log en tiempo real
tail -f logs/download_1kgp.log

# Ver archivos descargados
ls -lh data/vcf/ALL.chr*.vcf.gz | wc -l
```

---

## 🔧 PASO 2: CONSTRUIR PANELES (10 CORRIDAS CADA UNO)

Una vez completada la descarga:

```bash
./build_panels_master.sh
```

**Menú interactivo:**
```
1) Solo Panel A (3 superpoblaciones) - ~6-8 horas
2) Solo Panel B (N11, 11 poblaciones) - ~8-10 horas  
3) Ambos paneles (recomendado) - ~14-18 horas
```

### ¿Qué hace cada script?

#### `build_panel_3way_10runs.sh`

Para **cada corrida** (1-10):
1. Selecciona 100 muestras al azar de AFR, EUR, AMR (con reemplazo)
2. Para cada cromosoma (1-22):
   - Extrae esas muestras del VCF
   - Calcula frecuencias alélicas
3. Concatena los 22 cromosomas
4. Guarda como `PANEL_3WAY_run1.txt ... run10.txt`

#### `build_panel_N11_10runs.sh`

Para **cada corrida** (1-10):
1. Selecciona 100 muestras al azar de cada población (11 poblaciones)
2. Para cada cromosoma (1-22):
   - Extrae y calcula frecuencias
3. Concatena
4. Guarda como `PANEL_N11_run1.txt ... run10.txt`

---

## 📊 PASO 3: EJECUTAR iADMIX

### Opción A: Script automático (recomendado)

Crea el script ejecutor:

```bash
cat > run_iadmix_all.sh << 'EOF'
#!/bin/bash
set -euo pipefail

PROJ_DIR="$HOME/Documentos/Ancestria"
IADMIX="$PROJ_DIR/data/iadmix"
PANELS_DIR="$PROJ_DIR/refpanels"
RESULTS_DIR="$PROJ_DIR/results/iadmix_10runs"

mkdir -p "$RESULTS_DIR"

POOLS=("POOL1" "POOL2" "HOSPITAL")
BAMS=(
    "$PROJ_DIR/data/bam/01-50_samples_UDB-102_482263.merged.bam"
    "$PROJ_DIR/data/bam/51-100_samples_UDB-103_482264.merged.bam"
    "$PROJ_DIR/data/bam/HOSPITAL.bam"
)

echo "=========================================="
echo "EJECUTANDO iADMIX: 60 ANÁLISIS"
echo "=========================================="
echo ""

# Panel A: 3-way × 10 runs × 3 pools = 30 análisis
for run in {1..10}; do
    panel="$PANELS_DIR/PANEL_3WAY_run${run}.txt"
    
    for i in {0..2}; do
        pool="${POOLS[$i]}"
        bam="${BAMS[$i]}"
        output="$RESULTS_DIR/${pool}_3WAY_run${run}"
        
        echo "Panel 3WAY run$run - $pool"
        
        python2 "$IADMIX/runancestry.py" \
            -f "$panel" \
            --bam "$bam" \
            -p 100 \
            -o "$output" \
            --path "$IADMIX"
    done
done

# Panel B: N11 × 10 runs × 3 pools = 30 análisis
for run in {1..10}; do
    panel="$PANELS_DIR/PANEL_N11_run${run}.txt"
    
    for i in {0..2}; do
        pool="${POOLS[$i]}"
        bam="${BAMS[$i]}"
        output="$RESULTS_DIR/${pool}_N11_run${run}"
        
        echo "Panel N11 run$run - $pool"
        
        python2 "$IADMIX/runancestry.py" \
            -f "$panel" \
            --bam "$bam" \
            -p 100 \
            -o "$output" \
            --path "$IADMIX"
    done
done

echo ""
echo "✓ COMPLETADO: 60 análisis"
EOF

chmod +x run_iadmix_all.sh
./run_iadmix_all.sh
```

---

## 📈 PASO 4: ANALIZAR RESULTADOS

### Extraer resultados de las 10 corridas

```bash
cd ~/Documentos/Ancestria

# Panel 3WAY
for pool in POOL1 POOL2 HOSPITAL; do
    echo "=== $pool - Panel 3WAY ==="
    for run in {1..10}; do
        echo -n "Run$run: "
        grep "ADMIX_PROP" results/iadmix_10runs/${pool}_3WAY_run${run}.ancestry.out
    done
    echo ""
done

# Panel N11
for pool in POOL1 POOL2 HOSPITAL; do
    echo "=== $pool - Panel N11 ==="
    for run in {1..10}; do
        echo -n "Run$run: "
        grep "ADMIX_PROP" results/iadmix_10runs/${pool}_N11_run${run}.ancestry.out
    done
    echo ""
done
```

### Calcular medias ± desviación estándar

Crear script de análisis:

```python
#!/usr/bin/env python3
"""
Calcular media ± DE de 10 corridas
"""

import re
import numpy as np
from pathlib import Path

results_dir = Path("results/iadmix_10runs")

def parse_admix_prop(line):
    """Extrae proporciones de ancestría"""
    pops = {}
    parts = line.strip().split()
    for part in parts:
        if ':' in part:
            pop, val = part.split(':')
            pops[pop] = float(val)
    return pops

# Panel 3WAY
for pool in ["POOL1", "POOL2", "HOSPITAL"]:
    print(f"\n{pool} - Panel 3WAY:")
    
    all_runs = {"AFR": [], "EUR": [], "AMR": []}
    
    for run in range(1, 11):
        file = results_dir / f"{pool}_3WAY_run{run}.ancestry.out"
        
        with open(file) as f:
            for line in f:
                if "ADMIX_PROP" in line:
                    pops = parse_admix_prop(line)
                    for pop in ["AFR", "EUR", "AMR"]:
                        all_runs[pop].append(pops[pop])
    
    for pop in ["AFR", "EUR", "AMR"]:
        mean = np.mean(all_runs[pop])
        std = np.std(all_runs[pop])
        print(f"  {pop}: {mean:.4f} ± {std:.4f}")

# Repetir para N11...
```

---

## 📁 ESTRUCTURA FINAL

```
Ancestria/
├── data/
│   ├── vcf/
│   │   └── ALL.chr{1..22}.*.vcf.gz  (500GB)
│   └── panel/
│       └── integrated_call_samples_v3.20130502.ALL.panel
├── refpanels/
│   ├── PANEL_3WAY_run{1..10}.txt
│   └── PANEL_N11_run{1..10}.txt
└── results/
    └── iadmix_10runs/
        ├── POOL1_3WAY_run{1..10}.ancestry.out
        ├── POOL1_N11_run{1..10}.ancestry.out
        ├── POOL2_3WAY_run{1..10}.ancestry.out
        ├── POOL2_N11_run{1..10}.ancestry.out
        ├── HOSPITAL_3WAY_run{1..10}.ancestry.out
        └── HOSPITAL_N11_run{1..10}.ancestry.out
```

---

## ⚠️ SOLUCIÓN DE PROBLEMAS

### Error: "Espacio insuficiente"
```bash
# Limpiar archivos temporales
rm -rf work/panel_*_runs/

# Ver espacio usado
du -sh data/vcf
```

### Error: "No se encuentra bcftools"
```bash
# Instalar bcftools
conda install -c bioconda bcftools
# O
sudo apt install bcftools
```

### Descarga se interrumpe
El script tiene reintentos automáticos. Si se detiene completamente:
```bash
# Reanudar desde donde quedó
./download_1kgp_phase3.sh
```

---

## 📊 RESULTADOS ESPERADOS

### Panel 3WAY (AFR, EUR, AMR)
```
POOL1:     AFR: 52.8% ± 1.2%  EUR: 15.7% ± 0.8%  AMR: 31.5% ± 1.5%
POOL2:     AFR: 53.0% ± 1.1%  EUR: 15.6% ± 0.9%  AMR: 31.4% ± 1.4%
HOSPITAL:  AFR: 53.8% ± 1.3%  EUR: 15.1% ± 0.7%  AMR: 31.1% ± 1.6%
```

Las desviaciones estándar **deben ser pequeñas** (<2%) si el muestreo es robusto.

---

## 🔬 VALIDACIÓN CIENTÍFICA

**Criterios de calidad:**
1. ✅ 10 corridas independientes (muestreo aleatorio)
2. ✅ 100 individuos por población
3. ✅ Genoma completo (no solo chr22)
4. ✅ Muestreo con reemplazo (evita sesgo)
5. ✅ Media ± DE reportada

**Para publicación:**
- Reportar media ± DE de las 10 corridas
- Mencionar: "10 corridas independientes con selección aleatoria de 100 individuos con reemplazo"
- Comparar Panel 3WAY vs N11 para validar consistencia

---

## 📧 SOPORTE

Si encuentras errores o tienes dudas sobre el pipeline, revisa:
1. Los logs en `logs/`
2. El espacio disponible en disco
3. Que bcftools funcione correctamente
