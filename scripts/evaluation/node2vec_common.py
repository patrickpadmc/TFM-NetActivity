#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): Node2Vec implementado a mano.

No hay paquete 'node2vec' ni 'nodevectors' instalado en el venv (solo
gensim, networkx, numba, scipy). Se implementan caminatas aleatorias
sesgadas de 2do orden (parametros p, q, Grover & Leskovec 2016) sobre
la representacion CSR del grafo, aceleradas con numba (compilacion JIT,
paralelas via prange), y se entrena el skip-gram con gensim.Word2Vec
sobre las caminatas generadas.

Motivo de la implementacion propia (no precomputar tablas alias por
arista, como hace la implementacion "clasica" de node2vec): con grado
medio 2917.5 y maximo 18,707 (por la densidad de coexpressdb),
precomputar tablas alias para las ~220M aristas dirigidas del grafo
simetrizado es computacionalmente inviable (séria O(aristas x grado^2)
en el peor caso). En su lugar, los pesos de transicion se calculan al
vuelo, solo para los pasos de caminata realmente visitados.
"""
import numpy as np
import pandas as pd
from numba import njit, prange
from scipy.sparse import csr_matrix
from gensim.models import Word2Vec


def log(msg):
    print(f"[node2vec] {msg}", flush=True)


def build_csr(train_graph_path):
    df = pd.read_csv(train_graph_path, sep="\t", usecols=["node1", "node2"], dtype=str)
    nodes = sorted(set(df["node1"]) | set(df["node2"]))
    idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)

    rows = df["node1"].map(idx).to_numpy()
    cols = df["node2"].map(idx).to_numpy()
    all_rows = np.concatenate([rows, cols]).astype(np.int32)
    all_cols = np.concatenate([cols, rows]).astype(np.int32)
    data = np.ones(len(all_rows), dtype=np.int8)
    A = csr_matrix((data, (all_rows, all_cols)), shape=(n_nodes, n_nodes))
    A.data[:] = 1
    A.sum_duplicates()
    A.data[:] = 1
    A.sort_indices()  # imprescindible: la busqueda binaria en _has_neighbor asume orden

    indptr = A.indptr.astype(np.int64)
    indices = A.indices.astype(np.int32)
    return indptr, indices, nodes, idx


@njit(cache=True)
def _has_neighbor(indices, start, end, target):
    lo, hi = start, end
    while lo < hi:
        mid = (lo + hi) // 2
        v = indices[mid]
        if v == target:
            return True
        elif v < target:
            lo = mid + 1
        else:
            hi = mid
    return False


@njit(cache=True)
def _weighted_choice(weights, rand_val):
    total = 0.0
    for w in weights:
        total += w
    r = rand_val * total
    cum = 0.0
    for i in range(len(weights)):
        cum += weights[i]
        if r <= cum:
            return i
    return len(weights) - 1


@njit(cache=True, parallel=True)
def generate_walks(indptr, indices, start_nodes, walk_length, num_walks, p, q, seed):
    n_starts = len(start_nodes)
    total_walks = n_starts * num_walks
    walks = np.full((total_walks, walk_length), -1, dtype=np.int64)

    for wi in prange(total_walks):
        node_i = wi % n_starts
        s = start_nodes[node_i]
        np.random.seed(seed + wi)

        walks[wi, 0] = s
        deg_s = indptr[s + 1] - indptr[s]
        if deg_s == 0:
            continue
        first_next_idx = indptr[s] + np.random.randint(0, deg_s)
        cur = indices[first_next_idx]
        walks[wi, 1] = cur
        prev = s

        for step in range(2, walk_length):
            start_c = indptr[cur]
            end_c = indptr[cur + 1]
            deg_c = end_c - start_c
            if deg_c == 0:
                break
            weights = np.empty(deg_c, dtype=np.float64)
            start_p = indptr[prev]
            end_p = indptr[prev + 1]
            for j in range(deg_c):
                nxt = indices[start_c + j]
                if nxt == prev:
                    weights[j] = 1.0 / p
                elif _has_neighbor(indices, start_p, end_p, nxt):
                    weights[j] = 1.0
                else:
                    weights[j] = 1.0 / q
            r = np.random.random()
            choice = _weighted_choice(weights, r)
            nxt_node = indices[start_c + choice]
            walks[wi, step] = nxt_node
            prev = cur
            cur = nxt_node

    return walks


def walks_to_sentences(walks):
    sentences = []
    for row in walks:
        toks = [str(x) for x in row if x >= 0]
        if len(toks) >= 2:
            sentences.append(toks)
    return sentences


def train_node2vec(train_graph_path, p, q, walk_length, num_walks, dim=128,
                    window=10, negative=5, epochs=1, workers=8, seed=42,
                    indptr=None, indices=None, nodes=None):
    """
    Si se pasan indptr/indices/nodes precomputados, se reutilizan (evita
    reconstruir el CSR en cada combinacion de p,q durante el grid).
    Devuelve (embeddings_dict, timings_dict).
    """
    import time
    t0 = time.time()
    if indptr is None:
        indptr, indices, nodes, _ = build_csr(train_graph_path)
    n_nodes = len(nodes)
    t_build = time.time() - t0

    t1 = time.time()
    start_nodes = np.arange(n_nodes, dtype=np.int64)
    walks = generate_walks(indptr, indices, start_nodes, walk_length, num_walks,
                            np.float64(p), np.float64(q), seed)
    t_walks = time.time() - t1

    t2 = time.time()
    sentences = walks_to_sentences(walks)
    model = Word2Vec(sentences=sentences, vector_size=dim, window=window,
                      min_count=0, sg=1, negative=negative, workers=workers,
                      epochs=epochs, seed=seed)
    t_w2v = time.time() - t2

    embeddings = {}
    n_missing = 0
    for i, node_id in enumerate(nodes):
        key = str(i)
        if key in model.wv:
            embeddings[node_id] = model.wv[key]
        else:
            n_missing += 1

    timings = {
        "t_build_csr": t_build, "t_walks": t_walks, "t_word2vec": t_w2v,
        "n_walks": int(walks.shape[0]), "walk_length": walk_length,
        "n_missing_embeddings": n_missing,
    }
    return embeddings, timings, (indptr, indices, nodes)
