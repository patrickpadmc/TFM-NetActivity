#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): Spectral Embedding, dimension 128.

Intento sobre el grafo de entrenamiento combinado COMPLETO
(graph_edges_train_only.tsv: 110,324,395 aristas, 75,128 nodos, grado
medio 2917.5, grado maximo 18,707). Si no es computacionalmente viable
(memoria insuficiente o no convergencia dentro del limite de tiempo del
job), este script fallara de forma explicita y ese fallo, junto con el
uso de memoria/tiempo hasta el punto de falla, ES la documentacion del
limite computacional exigida por el protocolo -- en ese caso el
siguiente paso es construir el subgrafo comun (definido en Seccion 5)
y reintentar sobre el.

Metodo: los k=dim autovectores correspondientes a los autovalores mas
PEQUENOS del Laplaciano normalizado simetrico L = I - D^-1/2 A D^-1/2
(embedding espectral estandar / Laplacian Eigenmaps). Se usa
scipy.sparse.linalg.eigsh en modo shift-invert con un sigma pequeno
(1e-5, no exactamente 0) para evitar la singularidad exacta que
introduciria el autovalor 0 del Laplaciano (el grafo tiene 23
componentes conexas, Seccion 3, por lo que el autovalor 0 esta
repetido y la factorizacion LU en sigma=0 exacto seria singular).

Uso:
    python3 spectral_embedding_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --out data/processed/evaluation/graph_v4/embeddings/spectral_v4_full
"""
import argparse
import resource
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh


def log(msg):
    print(f"[spectral] {msg}", flush=True)


def peak_mem_gb():
    # ru_maxrss esta en KB en Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--sigma", type=float, default=1e-5)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--maxiter", type=int, default=5000)
    args = ap.parse_args()

    t0 = time.time()
    log("Cargando grafo de entrenamiento combinado...")
    df = pd.read_csv(args.train_graph, sep="\t", usecols=["node1", "node2"], dtype=str)
    log(f"{len(df):,} aristas cargadas")

    nodes = sorted(set(df["node1"]) | set(df["node2"]))
    idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)
    log(f"{n_nodes:,} nodos unicos")

    rows = df["node1"].map(idx).to_numpy()
    cols = df["node2"].map(idx).to_numpy()
    all_rows = np.concatenate([rows, cols])
    all_cols = np.concatenate([cols, rows])
    data = np.ones(len(all_rows), dtype=np.float64)
    A = csr_matrix((data, (all_rows, all_cols)), shape=(n_nodes, n_nodes))
    A.data[:] = 1.0
    A.sum_duplicates()
    A.data[:] = 1.0
    log(f"Adyacencia construida ({A.nnz:,} no-ceros, "
        f"densidad={A.nnz / (n_nodes ** 2):.5f}). "
        f"Mem pico: {peak_mem_gb():.1f} GB. t={time.time() - t0:.1f}s")

    log("Calculando Laplaciano normalizado simetrico...")
    L = laplacian(A, normed=True)
    L = L.tocsc()
    log(f"Laplaciano construido ({L.nnz:,} no-ceros). "
        f"Mem pico: {peak_mem_gb():.1f} GB. t={time.time() - t0:.1f}s")

    log(f"Calculando {args.dim} autovectores mas pequenos via eigsh "
        f"(shift-invert, sigma={args.sigma}, tol={args.tol}, maxiter={args.maxiter})...")
    try:
        vals, vecs = eigsh(L, k=args.dim, sigma=args.sigma, which="LM",
                            tol=args.tol, maxiter=args.maxiter)
    except Exception as e:
        log(f"FALLO: {type(e).__name__}: {e}")
        log(f"Mem pico en el momento del fallo: {peak_mem_gb():.1f} GB. "
            f"t={time.time() - t0:.1f}s")
        raise

    log(f"Autovectores calculados. Mem pico: {peak_mem_gb():.1f} GB. "
        f"t={time.time() - t0:.1f}s")

    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]

    np.savez(args.out, node_ids=np.array(nodes, dtype=object),
             embeddings=vecs, eigenvalues=vals)
    log(f"Guardado -> {args.out}.npz")
    log(f"Tiempo total: {time.time() - t0:.1f}s. Mem pico final: {peak_mem_gb():.1f} GB")


if __name__ == "__main__":
    main()
