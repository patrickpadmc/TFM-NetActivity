#!/usr/bin/env python
"""
Verifica el numero de componentes conexas del grafo unificado.

Diagnostico necesario tras observar un autovalor ~0 en el segundo
puesto del spectral embedding (k=256) - posible indicio de que el
grafo no es una sola componente conexa.

Uso:
    python check_connectivity.py --edges <path> --out <path>
"""
import argparse
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
import pyarrow.csv as pv


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", required=True)
    parser.add_argument("--out", required=True, help="Path al TSV de resumen de salida")
    args = parser.parse_args()

    log("Leyendo aristas...")
    convert_opts = pv.ConvertOptions(include_columns=["node1", "node2"])
    parse_opts = pv.ParseOptions(delimiter="\t")
    table = pv.read_csv(args.edges, parse_options=parse_opts, convert_options=convert_opts)
    node1 = table.column("node1").to_numpy(zero_copy_only=False)
    node2 = table.column("node2").to_numpy(zero_copy_only=False)
    log(f"Cargadas {len(node1)} aristas")

    log("Factorizando ids de nodo...")
    combined = np.concatenate([node1, node2])
    codes, uniques = pd.factorize(combined, sort=True)
    del combined
    n_edges = len(node1)
    row, col = codes[:n_edges], codes[n_edges:]
    n_nodes = len(uniques)
    log(f"{n_nodes} nodos unicos")

    log("Construyendo adyacencia sparse...")
    data = np.ones(len(row), dtype=np.float32)
    A = sp.coo_matrix((data, (row, col)), shape=(n_nodes, n_nodes)).tocsr()
    A = A.maximum(A.T)
    del row, col, data

    log("Calculando componentes conexas...")
    n_components, labels = connected_components(A, directed=False)
    log(f"Numero de componentes conexas: {n_components}")

    lines = [f"n_nodes\t{n_nodes}", f"n_components\t{n_components}"]

    if n_components > 1:
        sizes = np.bincount(labels)
        sizes_sorted = np.sort(sizes)[::-1]
        log(f"Tamano de las 10 componentes mas grandes: {sizes_sorted[:10]}")
        log(f"Nodos en la componente gigante: {sizes_sorted[0]} ({100*sizes_sorted[0]/n_nodes:.2f}%)")
        log(f"Nodos fuera de la componente gigante: {n_nodes - sizes_sorted[0]}")
        lines.append(f"giant_component_size\t{sizes_sorted[0]}")
        lines.append(f"giant_component_pct\t{100*sizes_sorted[0]/n_nodes:.4f}")
        lines.append("top10_component_sizes\t" + ",".join(map(str, sizes_sorted[:10])))

        labels_path = args.out.replace(".tsv", "_node_labels.tsv")
        pd.DataFrame({"node": uniques, "component": labels}).to_csv(labels_path, sep="\t", index=False)
        log(f"Etiquetas de componente por nodo guardadas en {labels_path}")
    else:
        lines.append(f"giant_component_size\t{n_nodes}")
        lines.append("giant_component_pct\t100.0")

    with open(args.out, "w") as f:
        f.write("\n".join(lines) + "\n")
    log(f"Resumen guardado en {args.out}")
    log("Listo")


if __name__ == "__main__":
    main()