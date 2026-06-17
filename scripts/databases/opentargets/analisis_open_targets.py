#!/usr/bin/python3

"""
Análisis de Open Targets Platform 26.03
Lee archivos Parquet y genera tabla de descriptivos target-disease-score
"""
 
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import pyarrow.parquet as pq
 
def analisis_open_targets(parquet_dir, output_dir="./opentargets_analisis"):
    """
    Analiza archivos Parquet de Open Targets y genera descriptivos
    
    Args:
        parquet_dir: directorio con archivos .parquet
        output_dir: directorio para guardar resultados
    """

    print("=" * 60)
    print("ANALISIS OPEN TARGETS PLATFORM 26.03")
    print("=" * 60)
    
    # Crear directorio de salida
    Path(output_dir).mkdir(exist_ok=True)
    
    # Encontrar y leer archivos Parquet
    print(f"\n[1/5] Buscando archivos Parquet en: {parquet_dir}")
    
    parquet_path = Path(parquet_dir)
    parquet_files = list(parquet_path.glob("**/*.parquet"))

    if not parquet_files:
        print(f" No se encontraron archivos .parquet en {parquet_dir}")
        return None
    
    print(f" Encontrados {len(parquet_files)} archivos parquet")
    
    # Cargar todos los parquets
    print("\n[2/5] Leyendo archivos Parquet...")
    
    dfs = []
    for i, pf in enumerate(parquet_files[:10]):  # Limitar a los primeros 10 para evitar memory issues
        try:
            table = pq.read_table(pf)
            df_temp = table.to_pandas()
            dfs.append(df_temp)
            print(f" {pf.name}: {len(df_temp):,} filas")
        except Exception as e:
            print(f" Error leyendo {pf.name}: {e}")
    
    if not dfs:
        print(" No se pudieron leer los archivos parquet")
        return None
    
    # Concatenar todos
    print(f"\n Combinando {len(dfs)} archivos...")
    df = pd.concat(dfs, ignore_index=True)
    print(f" Total: {len(df):,} registros")
    
    print(f" Columnas encontradas: {list(df.columns)}")
    
    
    # Identificar columnas clave
   
    print("\n[3/5] Identificando columnas clave...")
    
    # Buscar columnas de target, disease y score
    target_col = None
    disease_col = None
    score_col = None
    
    for col in df.columns:
        if 'target' in col.lower():
            target_col = col
            print(f" Target column: {target_col}")
        elif 'disease' in col.lower():
            disease_col = col
            print(f" Disease column: {disease_col}")
        elif 'score' in col.lower():
            score_col = col
            print(f" Score column: {score_col}")
    
    if not all([target_col, disease_col, score_col]):
        print("\n Algunas columnas clave no encontradas. Listando columnas disponibles:")
        for i, col in enumerate(df.columns, 1):
            print(f" {i}. {col}")
    
    
    # Estadisticas basicas
   
    print("\n[4/5] Calculando estadísticas básicas...")
    
    n_asociaciones = len(df)
    n_targets = df[target_col].nunique() if target_col else 0
    n_diseases = df[disease_col].nunique() if disease_col else 0
    
    # Score statistics
    score_stats = {}
    if score_col and df[score_col].dtype in [np.float64, np.float32]:
        score_stats = {
            'mean': df[score_col].mean(),
            'median': df[score_col].median(),
            'std': df[score_col].std(),
            'min': df[score_col].min(),
            'max': df[score_col].max(),
        }
    
    # Crear tabla de descriptivos
    descriptivos = {
        'Métrica': [
            'Total de asociaciones',
            'Targets únicos',
            'Enfermedades únicas',
            'Asociaciones por target (media)',
            'Asociaciones por disease (media)',
            'Score medio',
            'Score mediano',
            'Score std',
            'Score mín',
            'Score máx'
        ],
        'Valor': [
            f"{n_asociaciones:,}",
            f"{n_targets:,}",
            f"{n_diseases:,}",
            f"{n_asociaciones / max(n_targets, 1):.2f}",
            f"{n_asociaciones / max(n_diseases, 1):.2f}",
            f"{score_stats.get('mean', 'N/A'):.4f}" if score_stats else 'N/A',
            f"{score_stats.get('median', 'N/A'):.4f}" if score_stats else 'N/A',
            f"{score_stats.get('std', 'N/A'):.4f}" if score_stats else 'N/A',
            f"{score_stats.get('min', 'N/A'):.4f}" if score_stats else 'N/A',
            f"{score_stats.get('max', 'N/A'):.4f}" if score_stats else 'N/A',
        ]
    }
    
    descriptivos_df = pd.DataFrame(descriptivos)
    print("\n Descriptivos principales:")
    print(descriptivos_df.to_string(index=False))
    
    # Guardar descriptivos
    descriptivos_df.to_csv(f"{output_dir}/descriptivos.csv", index=False)
    print(f"\n Guardado en: {output_dir}/descriptivos.csv")
    

    # Tabla Target-Disease
    print("\n[5/5] Preparando tabla target-disease-score...")
    
    if target_col and disease_col and score_col:
        # Seleccionar columnas clave
        tds_tabla = df[[target_col, disease_col, score_col]].copy()
        tds_tabla.columns = ['target_id', 'disease_id', 'score']
        
        # Remover duplicados (si los hay)
        tds_tabla = tds_tabla.drop_duplicates()
        
        # Ordenar por score descendente
        tds_tabla = tds_tabla.sort_values('score', ascending=False)
        
        print(f" {len(tds_tabla):,} asociaciones únicas")
        print(f"\n Top 10 asociaciones (highest score):")
        print(tds_tabla.head(10).to_string(index=False))
        
        # Guardar tabla completa (o chunked si es muy grande)
        if len(tds_tabla) > 1000000:
            print(f"\n Tabla muy grande ({len(tds_tabla):,} filas)")
            print(f" Guardando top 500,000...")
            tds_tabla.head(500000).to_csv(f"{output_dir}/target_disease_score_top500k.csv", index=False)
            print(f" Guardado en: {output_dir}/target_disease_score_top500k.csv")
        else:
            tds_tabla.to_csv(f"{output_dir}/target_disease_score.csv", index=False)
            print(f" Tabla completa guardada en: {output_dir}/target_disease_score.csv")
        
        # Guardar top 1000
        tds_tabla.head(1000).to_csv(f"{output_dir}/target_disease_score_top1000.csv", index=False)
        print(f" Top 1000 guardado en: {output_dir}/target_disease_score_top1000.csv")
    

    # Distribucion de scores

    if score_col:
        print("\n Analizando distribución de scores...")
        bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        hist, bin_edges = np.histogram(df[score_col], bins=bins)
        
        distribution_data = {
            'Score Range': [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins)-1)],
            'Count': hist,
            'Percentage': [f"{(h/len(df)*100):.2f}%" for h in hist]
        }
        
        distribution_df = pd.DataFrame(distribution_data)
        print("\n Distribución de scores:")
        print(distribution_df.to_string(index=False))
        
        distribution_df.to_csv(f"{output_dir}/score_distribution.csv", index=False)
    
    # Resumen final
    print("\n" + "=" * 60)
    print("ANALISIS COMPLETADO")
    print("=" * 60)
    print(f"\nArchivos generados en: {output_dir}/")
    print("  1. descriptives.csv - Estadísticas principales")
    print("  2. target_disease_score.csv - Tabla completa (target-disease-score)")
    print("  3. target_disease_score_top1000.csv - Top 1000 asociaciones")
    print("  4. score_distribution.csv - Distribución de scores")
    print("\n Listo para análisis posterior\n")
    
    return tds_tabla if 'tds_tabla' in locals() else None
 
if __name__ == "__main__":
    if len(sys.argv) < 2:
        parquet_dir = "./association_by_datasource_direct"
        output_dir = "./opentargets_analysis"
        print(f"Uso: python3 analyze_opentargets.py <parquet_dir> [output_dir]")
        print(f"\nUsando valores por defecto:")
        print(f"  Input: {parquet_dir}")
        print(f"  Output: {output_dir}")
    else:
        parquet_dir = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "./opentargets_analisis"
    
    analisis_open_targets(parquet_dir, output_dir)
 