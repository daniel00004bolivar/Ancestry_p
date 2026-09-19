#!/usr/bin/env python3
"""
ANÁLISIS ESTADÍSTICO: 10 CORRIDAS DE ANCESTRÍA
Calcula media ± desviación estándar para cada pool y panel
"""

import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Configuración
RESULTS_DIR = Path("results/iadmix_10runs")
POOLS = ["POOL1", "POOL2", "HOSPITAL"]
PANELS = {
    "3WAY": ["AFR", "EUR", "AMR"],
    "N11": ["IBS", "CEU", "YRI", "ESN", "GWD", "LWK", "PEL", "MXL", "PUR", "FIN", "CLM"]
}

def parse_admix_line(line):
    """
    Extrae proporciones de ancestría de una línea ADMIX_PROP
    Ejemplo: "final maxval -123.45 ADMIX_PROP AFR:0.52 EUR:0.31 AMR:0.17"
    """
    pops = {}
    parts = line.strip().split()
    
    for part in parts:
        if ':' in part and not part.startswith('-'):
            try:
                pop, val = part.split(':')
                pops[pop] = float(val)
            except ValueError:
                continue
    
    return pops

def collect_results(pool, panel_name, n_runs=10):
    """
    Recolecta resultados de todas las corridas para un pool y panel
    """
    results = {pop: [] for pop in PANELS[panel_name]}
    
    for run in range(1, n_runs + 1):
        file_path = RESULTS_DIR / f"{pool}_{panel_name}_run{run}.ancestry.out"
        
        if not file_path.exists():
            print(f"⚠️  ADVERTENCIA: No se encuentra {file_path}")
            continue
        
        with open(file_path) as f:
            for line in f:
                if "ADMIX_PROP" in line:
                    pops = parse_admix_line(line)
                    for pop in PANELS[panel_name]:
                        if pop in pops:
                            results[pop].append(pops[pop])
    
    return results

def calculate_statistics(results):
    """
    Calcula media, desviación estándar, min, max
    """
    stats = {}
    
    for pop, values in results.items():
        if len(values) == 0:
            stats[pop] = {"mean": 0, "std": 0, "min": 0, "max": 0, "n": 0}
        else:
            stats[pop] = {
                "mean": np.mean(values),
                "std": np.std(values),
                "min": np.min(values),
                "max": np.max(values),
                "n": len(values)
            }
    
    return stats

def print_summary_table(pool, panel_name, stats):
    """
    Imprime tabla resumen bonita
    """
    print(f"\n{'='*70}")
    print(f"{pool} - Panel {panel_name}")
    print(f"{'='*70}")
    print(f"{'Población':<12} {'Media':<12} {'DE':<12} {'Min':<12} {'Max':<12} {'N'}")
    print('-'*70)
    
    for pop in PANELS[panel_name]:
        s = stats[pop]
        print(f"{pop:<12} {s['mean']:>10.4f}  {s['std']:>10.4f}  "
              f"{s['min']:>10.4f}  {s['max']:>10.4f}  {s['n']:>4}")

def export_to_csv(all_results, output_file):
    """
    Exporta resultados a CSV para análisis posterior
    """
    rows = []
    
    for (pool, panel), stats in all_results.items():
        for pop, s in stats.items():
            rows.append({
                "Pool": pool,
                "Panel": panel,
                "Población": pop,
                "Media": s["mean"],
                "DE": s["std"],
                "Min": s["min"],
                "Max": s["max"],
                "N_corridas": s["n"]
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    print(f"\n✓ Resultados exportados a: {output_file}")

def main():
    if not RESULTS_DIR.exists():
        print(f"ERROR: No se encuentra el directorio {RESULTS_DIR}")
        print("Ejecuta primero run_iadmix_all.sh")
        sys.exit(1)
    
    print("\n" + "="*70)
    print("ANÁLISIS ESTADÍSTICO: 10 CORRIDAS DE ANCESTRÍA")
    print("="*70)
    
    all_results = {}
    
    # Procesar cada combinación de pool y panel
    for pool in POOLS:
        for panel_name in PANELS.keys():
            print(f"\nProcesando {pool} - {panel_name}...")
            
            results = collect_results(pool, panel_name)
            stats = calculate_statistics(results)
            
            all_results[(pool, panel_name)] = stats
            print_summary_table(pool, panel_name, stats)
    
    # Exportar a CSV
    output_csv = "results/ancestry_statistics_10runs.csv"
    export_to_csv(all_results, output_csv)
    
    # Resumen final
    print("\n" + "="*70)
    print("RESUMEN PARA REPORTE")
    print("="*70)
    
    for pool in POOLS:
        print(f"\n{pool}:")
        
        # Panel 3WAY
        stats_3way = all_results[(pool, "3WAY")]
        afr = stats_3way["AFR"]["mean"] * 100
        afr_std = stats_3way["AFR"]["std"] * 100
        eur = stats_3way["EUR"]["mean"] * 100
        eur_std = stats_3way["EUR"]["std"] * 100
        amr = stats_3way["AMR"]["mean"] * 100
        amr_std = stats_3way["AMR"]["std"] * 100
        
        print(f"  Panel 3WAY:")
        print(f"    AFR: {afr:.2f}% ± {afr_std:.2f}%")
        print(f"    EUR: {eur:.2f}% ± {eur_std:.2f}%")
        print(f"    AMR: {amr:.2f}% ± {amr_std:.2f}%")
        
        # Panel N11 - Agregado a superpoblaciones
        stats_n11 = all_results[(pool, "N11")]
        
        afr_n11 = (stats_n11["YRI"]["mean"] + stats_n11["ESN"]["mean"] + 
                   stats_n11["GWD"]["mean"] + stats_n11["LWK"]["mean"]) * 100
        eur_n11 = (stats_n11["IBS"]["mean"] + stats_n11["CEU"]["mean"] + 
                   stats_n11["FIN"]["mean"]) * 100
        amr_n11 = (stats_n11["CLM"]["mean"] + stats_n11["PEL"]["mean"] + 
                   stats_n11["MXL"]["mean"] + stats_n11["PUR"]["mean"]) * 100
        
        print(f"  Panel N11 (agregado):")
        print(f"    AFR: {afr_n11:.2f}%")
        print(f"    EUR: {eur_n11:.2f}%")
        print(f"    AMR: {amr_n11:.2f}%")
        print(f"    (CLM específico: {stats_n11['CLM']['mean']*100:.2f}% ± "
              f"{stats_n11['CLM']['std']*100:.2f}%)")

if __name__ == "__main__":
    main()
