#!/usr/bin/python3

"""
Análisis de DisGeNET 25.4
Genera tabla de descriptivos y extrae gene-disease scores
"""

import pandas as pd
import numpy as np
import gzip
import sys
from pathlib import Path

def analisis_disgenet(input_file, output_dir="./disgenet_analisis"):
    """
    Analiza el archivo de DisGeNET y genera descriptivos

    Argumentos:
    input_file: ruta al archivo .tsv.gz de DisGeNET
    output_dir: directorio para guardar resultados
    """

    print("="*60)
    print('ANÁLISIS DisGeNET 25.4')
    print("="*60)

    # Crear directorio de salida
    Path(output_dir).mkdir(exist_ok=True)

    # Detectar si está comprimido
    if input_file.endswith('.gz'):
        print(f"n[1/5] Leyendo archivo comprimido: {input_file}")
        df = pd.read_csv(input_file, sep='\t', compression='gzip', low_memory=False)
    else:
        print(f"[1/5] Leyendo archivo: {input_file}")
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    print(f' {len(df):,} asociaciones gene-disease cargadas.')
    print(f' Columnas: {list(df.columns)}')


    # ESTADISITICAS BASISCAS
    print("\n[2/5] Calculando estadisticas basicas...")

    n_asociaciones = len(df)

    # Genes unicos
    if 'geneId' in df.columns or 'gene_id' in df.columns:
        gene_col = [c for c in df.columns if 'gene' in c.lower()][0]
        n_genes = df[gene_col].nunique()
    else:
        gene_col = None
        n_genes = 0

    # Enfermedades unicas
    if 'diseaseId' in df. columns or 'disease_id' in df.columns:
        disease_col = [c for c in df.columns if 'disease' in c.lower()][0]
        n_diseases = df[disease_col].nunique()
    else:
        disease_col = None
        n_diseases = 0

    # Score statitistics
    score_cols = [c for c in df.columns if 'score' in c.lower()]
    if score_cols:
        score_col = score_cols[0]
        score_stats = {
            'mean': df[score_col].mean(),
            'median': df[score_col].median(),
            'std': df[score_col].std(),
            'min': df[score_col].min(),
            'max': df[score_col].max(),
        }
    else:
        score_col = None
        score_stats = {}

        descriptivos = {
            'Métrica': [
                'Total de asociaciones',
                'Genes únicos',
                'Enfermedades únicas',
                'Associations por gene (media)',
                'Associations por disease (media)',
                'Score medio',
                'Score mediano',
                'Score std',
                'Score mín',
                'Score máx'
            ],
            'Valor': [
                f"{n_asociaciones:,}",
                f"{n_genes:,}",
                f"{n_diseases:,}",
                f"{n_asociaciones / max(n_genes, 1):.2f}",
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
        print(descriptivos_df.tostring(index=False))

        # Guardar descriptivos
        descriptivos_df.to_csv(f"{output_dir}/descriptivos.csv", index=False)
        print(f"\n Guardado en: {output_dir}/descriptivos.csv")

        # GENE-DISEASE SCORES
        print("\n[3/5] Preparando tabla gene-disease score...")

        if gene_col and disease_col and score_col:
            # Seleccionar columnas clave
            gds_tabla = df[[gene_col, disease_col, score_col]].copy()
            gds_tabla.columns = ['gene_id', 'disease_id', 'score']

            gds_tabla = gds_tabla.sort_values('score', ascending=False)

            print(f" {len(gds_tabla):,} asociaciones")
            print(f"\n Top 10 asosciaciones (highest scores):")
            print(gds_tabla.head(10).to_string(index=False))

            # Guardar tabla completa
            gds_tabla.to_csv(f"{output_dir}/gene_disease_scores.csv", index=False)
            print(f"\n Tabla completa guardada en:: {output_dir}/gene_disease_scores.csv")

            gds_tabla.head(1000).to_csv(f"{output_dir}/gene_disease_scores_top1000.csv", index=False)

            print(f"\n Top 1000 asociaciones guardadas en: {output_dir}/gene_disease_scores_top1000.csv")
        
        # Distribuucion de scores
        print("\n[4/5] Analizando distribución de scores...")

        if score_col:
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            hist, bin_edges = np.histogram(df[score_col], bins=bins)

            distribucion = {
                'Score Range': [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins)-1)],
                'Count':hist,
                'Percentage': [f"{(h/len(df)*100):.2f}%" for h in hist]
            }

            distribucion_df = pd.DataFrame(distribucion)
            print("\n Distribución de scores:")
            print(distribucion_df.to_string(index=False))

            distribucion_df.to_csv(f"{output_dir}/score_distribution.csv", index=False)

        # Analisis de fuentes 
        print("\n[5/5] Analizando por fuentes de datos...")

        source_cols = [c for c in df.columns if 'source' in c.lower()]
        if source_cols:
            source_col = source_cols[0]
            source_counts = df[source_col].value_counts()

            print("\n Asociaciones por fuente:")
            for source, count in source_counts.head(10).items():
                pct = (count/len(df))*100
                print(f" - {source}: {count:,} ({pct:.2f}%)")

            source_counts.to_frame('count').to_csv(f"{output_dir}/sources_distribution.csv")

    # Resumen final
    print("\n" + "=" * 60)
    print("ANÁLISIS COMPLETADO")
    print("=" * 60)
    print(f"\nArchivos generados en: {output_dir}/")
    print("  1. descriptives.csv - Estadísticas principales")
    print("  2. gene_disease_score.csv - Tabla completa (gene-disease-score)")
    print("  3. gene_disease_score_top1000.csv - Top 1000 asociaciones")
    print("  4. score_distribution.csv - Distribución de scores")
    print("  5. sources_distribution.csv - Asociaciones por fuente")
    print("\n Listo para análisis posterior\n")

    return gds_tabla if 'gds_tabla' in locals() else None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Valores por defecto (ajusta según tu estructura)
        input_file = "all_gene_disease_associations.tsv.gz"
        output_dir = "./disgenet_analisis"
        print(f"Uso: python3 analisis_disgenet.py <input_file> [output_dir]")
        print(f"\nUsando valores por defecto:")
        print(f"  Input: {input_file}")
        print(f"  Output: {output_dir}")
    else:
        input_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "./disgenet_analisis"
    
    analisis_disgenet(input_file, output_dir)