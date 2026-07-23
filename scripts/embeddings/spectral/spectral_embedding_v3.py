#!/usr/bin/env python
"""
Spectral embedding para graph_v3 de TFM-NetActivity.

Mismo approach que v2: Laplaciano normalizado del grafo (no dirigido,
no ponderado), autovectores de menor autovalor via LAPACK denso
(eigh con subset_by_index). Se descarta ARPACK sparse shift-invert por
la razon ya validada en v2: el Laplaciano tiene autovalor exacto en 0,
lo que rompe la factorizacion shift-invert cuando sigma=0.

Cambios respecto a v2:
  - Lee node1/node2 desde graph_v3/graph_edges.tsv (mismo esquema de
    columnas: node1, node2, node1_type, node2_type, relation, database;
    solo se usan node1/node2 para el spectral, igual que v2).
  - Nombres de salida versionados: spectral_embeddings_v3_k{k}.npy,
    spectral_eigenvalues_v3_k{k}.npy, node_order_v3.tsv.

Uso:
    python spectral_embedding_v3.py --edges <path> --out-dir <path> --k 256
"""
import argparse
import os
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.linalg as sla
import pyarrow.csv as pv


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_edges(edges_path):
    log(f"Leyendo aristas de {edges_path}")
    convert_opts = pv.ConvertOptions(include_columns=["node1", "node2"])
    parse_opts = pv.ParseOptions(delimiter="\t")
    table = pv.read_csv(edges_path, parse_options=parse_opts, convert_options=convert_opts)
    node1 = table.column("node1").to_numpy(zero_copy_only=False)
    node2 = table.column("node2").to_numpy(zero_copy_only=False)
    log(f"Cargadas {len(node1)} aristas")
    return node1, node2


def build_adjacency(node1, node2, out_dir):
    log("Factorizando ids de nodo")
    combined = np.concatenate([node1, node2])
    codes, uniques = pd.factorize(combined, sort=True)
    n_edges = len(node1)
    row = codes[:n_edges]
    col = codes[n_edges:]
    n_nodes = len(uniques)
    log(f"{n_nodes} nodos unicos")

    node_order_path = os.path.join(out_dir, "node_order_v3.tsv")
    pd.Series(uniques).to_csv(node_order_path, sep="\t", index=False, header=False)
    log(f"Guardado orden de nodos en {node_order_path} - fila i del embedding = este nodo")

    log("Construyendo matriz de adyacencia sparse")
    data = np.ones(len(row), dtype=np.float32)
    A = sp.coo_matrix((data, (row, col)), shape=(n_nodes, n_nodes)).tocsr()
    A.data[:] = 1.0
    A = A.maximum(A.T)
    A.setdiag(0)
    A.eliminate_zeros()
    density = A.nnz / (n_nodes ** 2)
    log(f"Adyacencia: shape={A.shape}, nnz={A.nnz}, densidad={density:.4f}")
    return A, n_nodes


def normalized_laplacian(A):
    log("Calculando Laplaciano normalizado")
    deg = np.asarray(A.sum(axis=1)).flatten()
    with np.errstate(divide="ignore"):
        deg_inv_sqrt = np.power(deg, -0.5)
    deg_inv_sqrt[np.isinf(deg_inv_sqrt)] = 0.0
    D_inv_sqrt = sp.diags(deg_inv_sqrt)
    n = A.shape[0]
    L = sp.eye(n, format="csr") - D_inv_sqrt @ A @ D_inv_sqrt
    return L


def eigh_dense(L, k_total):
    log("Convirtiendo Laplaciano sparse a denso")
    t0 = time.time()
    L_dense = np.asarray(L.todense(), dtype=np.float32)
    del L
    log(f"Conversion a denso completada en {time.time() - t0:.1f}s, shape={L_dense.shape}")

    log(f"Ejecutando LAPACK eigh (subset_by_index 0..{k_total - 1})")
    t0 = time.time()
    eigenvalues, eigenvectors = sla.eigh(L_dense, subset_by_index=[0, k_total - 1])
    log(f"Eigh denso completado en {time.time() - t0:.1f}s")
    return eigenvalues, eigenvectors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", required=True, help="Path a graph_edges.tsv de graph_v3")
    parser.add_argument("--out-dir", required=True, help="Directorio de salida")
    parser.add_argument("--k", type=int, default=256, help="Dimension del embedding")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    node1, node2 = load_edges(args.edges)
    A, n_nodes = build_adjacency(node1, node2, args.out_dir)
    L = normalized_laplacian(A)
    del A

    k = args.k
    k_total = k + 1  # incluye el autovector trivial (autovalor ~0), se descarta luego

    eigenvalues, eigenvectors = eigh_dense(L, k_total)

    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    embeddings = eigenvectors[:, 1:k_total]
    used_eigenvalues = eigenvalues[1:k_total]

    emb_path = os.path.join(args.out_dir, f"spectral_embeddings_v3_k{k}.npy")
    eig_path = os.path.join(args.out_dir, f"spectral_eigenvalues_v3_k{k}.npy")
    np.save(emb_path, embeddings)
    np.save(eig_path, used_eigenvalues)

    log(f"Embeddings guardados {embeddings.shape} en {emb_path}")
    log(f"Autovalores guardados {used_eigenvalues.shape} en {eig_path}")
    log("Listo")


if __name__ == "__main__":
    main()
