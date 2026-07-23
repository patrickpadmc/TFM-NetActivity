#!/usr/bin/env python
"""
Construye el archivo .edg (tab-separated, sin header) que pecanpy necesita
para node2vec, a partir de graph_v3/graph_edges.tsv.

El edge file de graph_v3 tiene aristas repetidas entre el mismo par de
nodos (misma relacion en distintas databases, o distintas relaciones
entre el mismo par). Para node2vec (grafo no dirigido, no ponderado)
esto es redundante: cada par (a,b) debe aparecer una sola vez. Este
script deduplica canonicalizando cada par sin importar el orden
(a,b) == (b,a), usando factorize + encoding a un solo entero para que
el unique sea rapido incluso con 110M filas de entrada.

Uso:
    python build_edgelist_v3.py --edges <path a graph_edges.tsv> --out <path al .edg de salida>
"""
import argparse
import time
import numpy as np
import pyarrow.csv as pv


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", required=True, help="Path a graph_edges.tsv de graph_v3")
    parser.add_argument("--out", required=True, help="Path de salida del .edg deduplicado")
    args = parser.parse_args()

    log(f"Leyendo aristas de {args.edges}")
    convert_opts = pv.ConvertOptions(include_columns=["node1", "node2"])
    parse_opts = pv.ParseOptions(delimiter="\t")
    table = pv.read_csv(args.edges, parse_options=parse_opts, convert_options=convert_opts)
    node1 = table.column("node1").to_numpy(zero_copy_only=False)
    node2 = table.column("node2").to_numpy(zero_copy_only=False)
    n_edges = len(node1)
    log(f"Cargadas {n_edges} filas")

    log("Factorizando ids de nodo")
    import pandas as pd
    combined = np.concatenate([node1, node2])
    codes, uniques = pd.factorize(combined, sort=True)
    n_nodes = len(uniques)
    row = codes[:n_edges].astype(np.int64)
    col = codes[n_edges:].astype(np.int64)
    log(f"{n_nodes} nodos unicos")

    log("Canonicalizando pares no dirigidos y removiendo self-loops")
    lo = np.minimum(row, col)
    hi = np.maximum(row, col)
    mask = lo != hi
    lo, hi = lo[mask], hi[mask]
    log(f"Removidos {(~mask).sum()} self-loops")

    log("Deduplicando pares (encoding a entero unico)")
    key = lo * n_nodes + hi
    unique_key = np.unique(key)
    lo_u = unique_key // n_nodes
    hi_u = unique_key % n_nodes
    log(f"Pares unicos: {len(unique_key)} (de {mask.sum()} tras remover self-loops)")

    log(f"Escribiendo edgelist en {args.out}")
    out_df = pd.DataFrame({"node1": uniques[lo_u], "node2": uniques[hi_u]})
    out_df.to_csv(args.out, sep="\t", header=False, index=False)

    log("Listo")


if __name__ == "__main__":
    main()
