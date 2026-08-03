import pandas as pd
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 20)
df = pd.read_csv('data/processed/analysis/seccion11_comparacion_embeddings/seccion11_metricas_comparacion.tsv', sep='\t')
cols = ['modelo','externo','n_genes_comunes','spearman','jaccard_k10','jaccard_k25','jaccard_k50','cka_lineal','ridge_r2']
print(df[cols].to_string(index=False))
