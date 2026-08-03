#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): baseline no aprendido.

Tres heuristicas simples, sin ningun parametro entrenado, calculadas
EXCLUSIVAMENTE sobre graph_edges_train_only.tsv (nunca sobre val/test,
para no filtrar informacion de evaluacion dentro del propio baseline):

  - Preferential attachment (PA):  score(u,v) = deg(u) * deg(v)
  - Common neighbors (CN):         score(u,v) = |N(u) intersec N(v)|
  - Jaccard:                       score(u,v) = CN(u,v) / |N(u) union N(v)|

El grafo se trata como no dirigido y agnostico al tipo de nodo para el
calculo de vecinos (N(u) incluye vecinos de cualquier tipo, alcanzados
por cualquier relacion) -- esto cubre de forma uniforme tanto la unica
relacion homogenea (GEN-ppi-GEN, donde "common neighbors" es la
definicion de manual) como las 6 relaciones bipartitas (donde Jaccard
sobre vecinos globales es la forma estandar de ponderar una proyeccion
bipartita simple, sin necesitar una implementacion aparte).

Al no haber entrenamiento, el score crudo de cada heuristica se evalua
directamente contra el test de cada tarea (AUPRC/AUROC via
sklearn.metrics), sin necesidad de la particion de validacion.

Uso:
    python3 run_baseline_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out data/processed/evaluation/graph_v4/baseline_results_v4.tsv
"""
import argparse
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.metrics import average_precision_score, roc_auc_score

TASKS = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
         "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]


def log(msg):
    print(f"[baseline] {msg}", flush=True)


def build_adjacency(train_graph_path):
    log("Cargando grafo de entrenamiento combinado...")
    df = pd.read_csv(train_graph_path, sep="\t", usecols=["node1", "node2"], dtype=str)
    log(f"{len(df):,} aristas cargadas")

    nodes = sorted(set(df["node1"]) | set(df["node2"]))
    idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)
    log(f"{n_nodes:,} nodos unicos (indice construido)")

    rows = df["node1"].map(idx).to_numpy()
    cols = df["node2"].map(idx).to_numpy()
    # Simetrizamos: grafo no dirigido para el calculo de vecinos/grado.
    all_rows = np.concatenate([rows, cols])
    all_cols = np.concatenate([cols, rows])
    data = np.ones(len(all_rows), dtype=np.int8)
    A = csr_matrix((data, (all_rows, all_cols)), shape=(n_nodes, n_nodes))
    A.data[:] = 1  # colapsa duplicados (p.ej. relaciones repetidas entre el mismo par) a presencia binaria
    A.sum_duplicates()
    A.data[:] = 1

    deg = np.asarray(A.sum(axis=1)).ravel()
    log(f"Grado minimo={deg.min()}, maximo={deg.max()}, medio={deg.mean():.1f}")
    return A, idx, deg


def score_pairs(A, deg, idx, pairs_df, batch_size=5000):
    """
    Calcula PA/CN/Jaccard en lotes. Necesario porque indexar la matriz
    dispersa con las filas de TODOS los pares de una tarea de una sola
    vez (A[us], con us de cientos de miles de elementos y filas de
    hasta 18,707 no-ceros repetidas muchas veces) puede generar una
    matriz intermedia con mas de 2^31 no-ceros, superando el limite de
    un entero de 32 bits que usa scipy internamente para los indices
    de la matriz dispersa -- eso corrompe memoria (segfault observado
    en el primer intento, job 308472, con batch_size implicito = todo
    el test set de una vez). Procesando en lotes acotados se evita.
    """
    valid = pairs_df["n1"].isin(idx) & pairs_df["n2"].isin(idx)
    n_dropped = int((~valid).sum())
    if n_dropped:
        log(f"  AVISO: {n_dropped} pares con algun nodo ausente del grafo de entrenamiento, descartados")
    pairs_df = pairs_df[valid]

    us_all = pairs_df["n1"].map(idx).to_numpy()
    vs_all = pairs_df["n2"].map(idx).to_numpy()
    n = len(us_all)

    pa = np.empty(n, dtype=np.float64)
    cn = np.empty(n, dtype=np.float64)
    n_batches = (n + batch_size - 1) // batch_size
    for bi, start in enumerate(range(0, n, batch_size)):
        end = min(start + batch_size, n)
        us = us_all[start:end]
        vs = vs_all[start:end]
        pa[start:end] = deg[us].astype(np.float64) * deg[vs].astype(np.float64)
        Au = A[us]
        Av = A[vs]
        cn[start:end] = np.asarray(Au.multiply(Av).sum(axis=1)).ravel().astype(np.float64)
        if (bi + 1) % 20 == 0 or (bi + 1) == n_batches:
            log(f"    lote {bi + 1}/{n_batches}")

    union = deg[us_all].astype(np.float64) + deg[vs_all].astype(np.float64) - cn
    jaccard = np.divide(cn, union, out=np.zeros_like(cn), where=union > 0)

    return pa, cn, jaccard, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-graph", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    t0 = time.time()
    A, idx, deg = build_adjacency(args.train_graph)
    log(f"Adyacencia construida en {time.time() - t0:.1f}s")

    results = []
    for task in TASKS:
        for split_type, test_name, neg_name in [
            ("edge", "test.tsv", "negatives_test.tsv"),
        ]:
            t1 = time.time()
            task_dir = f"{args.eval_dir}/{task}"
            pos = pd.read_csv(f"{task_dir}/{test_name}", sep="\t", dtype=str)
            neg = pd.read_csv(f"{task_dir}/{neg_name}", sep="\t", dtype=str)
            pairs = pd.concat([pos, neg], ignore_index=True)
            y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])

            pa, cn, jac, n_used = score_pairs(A, deg, idx, pairs)
            valid_mask = (pairs["n1"].isin(idx) & pairs["n2"].isin(idx)).to_numpy()
            y_aligned = y[valid_mask]

            for method_name, scores in [("baseline_pa", pa), ("baseline_cn", cn), ("baseline_jaccard", jac)]:
                auprc = average_precision_score(y_aligned, scores)
                auroc = roc_auc_score(y_aligned, scores)
                results.append({
                    "task": task, "split_type": split_type, "method": method_name,
                    "auprc": auprc, "auroc": auroc,
                    "n_test": int(len(y_aligned)), "n_pos_test": int(y_aligned.sum()),
                })
                log(f"  {task}/{split_type}/{method_name}: AUPRC={auprc:.4f} AUROC={auroc:.4f} "
                    f"(n_test={len(y_aligned):,}, n_pos={int(y_aligned.sum()):,})")
            log(f"  {task}/{split_type} evaluado en {time.time() - t1:.1f}s")

    # Split estricto (cold-start) de GEN-ass-DIS
    t1 = time.time()
    strict_dir = f"{args.eval_dir}/GEN-ass-DIS"
    pos = pd.read_csv(f"{strict_dir}/test_strict.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{strict_dir}/negatives_test_strict.tsv", sep="\t", dtype=str)
    pairs = pd.concat([pos, neg], ignore_index=True)
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    pa, cn, jac, n_used = score_pairs(A, deg, idx, pairs)
    valid_mask = (pairs["n1"].isin(idx) & pairs["n2"].isin(idx)).to_numpy()
    y_aligned = y[valid_mask]
    for method_name, scores in [("baseline_pa", pa), ("baseline_cn", cn), ("baseline_jaccard", jac)]:
        auprc = average_precision_score(y_aligned, scores)
        auroc = roc_auc_score(y_aligned, scores)
        results.append({
            "task": "GEN-ass-DIS", "split_type": "node_strict", "method": method_name,
            "auprc": auprc, "auroc": auroc,
            "n_test": int(len(y_aligned)), "n_pos_test": int(y_aligned.sum()),
        })
        log(f"  GEN-ass-DIS/node_strict/{method_name}: AUPRC={auprc:.4f} AUROC={auroc:.4f} "
            f"(n_test={len(y_aligned):,}, n_pos={int(y_aligned.sum()):,})")
    log(f"  GEN-ass-DIS/node_strict evaluado en {time.time() - t1:.1f}s")

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, sep="\t", index=False)
    log(f"Resultados -> {args.out}")
    log(f"Tiempo total: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
