#!/usr/bin/env python3
"""Consolida las 11 tablas __svm_results.tsv de Seccion 12 en una sola
tabla embedding_source x relation x metric, y genera la matriz de AUPRC
(fuente x relacion) mas un resumen de hiperparametros/umbrales."""
import glob
import pandas as pd

OUT_DIR = "data/processed/analysis/seccion12_svm_link_prediction"

files = sorted(glob.glob(f"{OUT_DIR}/*__svm_results.tsv"))
print(f"Archivos encontrados: {len(files)}")
dfs = [pd.read_csv(f, sep="\t") for f in files]
full = pd.concat(dfs, ignore_index=True)
print(f"Filas totales: {len(full)} (esperado: 9 fuentes x 7 relaciones + 2 fuentes x 1 relacion = 65)")

full_path = f"{OUT_DIR}/seccion12_tabla_completa.tsv"
full.to_csv(full_path, sep="\t", index=False)
print(f"Tabla completa -> {full_path}")

pivot_auprc = full.pivot(index="source", columns="relation", values="test_auprc")
pivot_path = f"{OUT_DIR}/seccion12_matriz_auprc.tsv"
pivot_auprc.to_csv(pivot_path, sep="\t")
print(f"Matriz AUPRC -> {pivot_path}")

hp_cols = ["source", "relation", "modelo_elegido", "linear_best_c", "linear_val_auprc",
           "rbf_best_c", "rbf_best_gamma", "rbf_val_auprc", "threshold", "f1_val_en_umbral", "class_weight"]
hp_path = f"{OUT_DIR}/seccion12_hiperparametros_umbrales.tsv"
full[hp_cols].to_csv(hp_path, sep="\t", index=False)
print(f"Hiperparametros/umbrales -> {hp_path}")

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
print("\n--- matriz AUPRC (fuente x relacion) ---")
print(pivot_auprc.round(3).to_string())

print("\n--- resumen por fuente (media AUPRC/AUROC/Hits@10/MRR sobre relaciones evaluadas) ---")
summary = full.groupby("source")[["test_auprc", "test_auroc", "hits_at_10", "mrr"]].mean().round(3)
print(summary.to_string())

print("\nListo.")
