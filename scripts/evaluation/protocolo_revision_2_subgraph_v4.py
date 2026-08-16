#!/usr/bin/env python3
import json

p = "data/processed/evaluation/graph_v4/common_subgraph_report.json"
with open(p) as f:
    r = json.load(f)

for k in [
    "seed", "metodo_muestreo", "fracciones",
    "n_aristas_grafo_entrenamiento", "n_aristas_subgrafo",
    "fraccion_aristas_global", "n_nodos_grafo_v4",
    "n_nodos_subgrafo", "n_aristas_protegidas_grado1",
    "n_nodos_rescatados"
]:
    print(f"{k}={r[k]}")

print("\nCobertura por tipo:")
for t, x in r["cobertura_por_tipo_nodo"].items():
    print(f"{t}: {x}")

print("\nRelaciones primarias:")
for x in r["detalle_por_relacion"]:
    if x["categoria"] == "primaria":
        print(x)
