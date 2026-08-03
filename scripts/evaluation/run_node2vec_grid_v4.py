#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): busqueda de p,q para Node2Vec.

Evalua las 16 combinaciones p,q en {0.25, 0.5, 1, 2} x {0.25, 0.5, 1, 2},
SOLO con AUPRC de validacion (metrica primaria), usando un presupuesto
de caminatas REDUCIDO (num_walks=3, walk_length=40 -- decision de costo
computacional confirmada: ~7x mas barato que el presupuesto completo
num_walks=10/walk_length=80 que se reserva para la configuracion final
evaluada en test). El presupuesto reducido se usa unicamente para
RANKEAR p,q; no se reporta como resultado final del TFM.

La seleccion de mejor (p,q) se hace POR TAREA (no una unica eleccion
global), ya que el mismo embedding Node2Vec se evalua de forma
independiente en cada una de las 7 tareas primarias, y nada garantiza
que el mismo (p,q) sea optimo para todas.

Uso:
    python3 run_node2vec_grid_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/evaluation/graph_v4/node2vec_grid
"""
import argparse
import itertools
import json
import time

import numpy as np
import pandas as pd

from node2vec_common import build_csr, train_node2vec
from common_eval import build_features, train_classifier

TASKS = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
         "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]

P_GRID = [0.25, 0.5, 1, 2]
Q_GRID = [0.25, 0.5, 1, 2]

GRID_WALK_LENGTH = 40
GRID_NUM_WALKS = 3
DIM = 128
SEED = 42


def log(msg):
    print(f"[grid] {msg}", flush=True)


def load_pairs_labels(task_dir, split):
    pos = pd.read_csv(f"{task_dir}/{split}.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}.tsv", sep="\t", dtype=str)
    pairs = list(zip(pos["n1"], pos["n2"])) + list(zip(neg["n1"], neg["n2"]))
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    return pairs, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-graph", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    t0 = time.time()
    log("Construyendo CSR del grafo de entrenamiento (reutilizado en las 16 combinaciones)...")
    indptr, indices, nodes, _ = build_csr(args.train_graph)
    log(f"{len(nodes):,} nodos. CSR listo en {time.time() - t0:.1f}s")

    # Cache de pares/labels por tarea (train y val), reutilizado en las 16 combinaciones.
    task_data = {}
    for task in TASKS:
        task_dir = f"{args.eval_dir}/{task}"
        train_pairs, train_y = load_pairs_labels(task_dir, "train")
        val_pairs, val_y = load_pairs_labels(task_dir, "val")
        task_data[task] = (train_pairs, train_y, val_pairs, val_y)
        log(f"  {task}: train={len(train_pairs):,} val={len(val_pairs):,}")

    results = []
    for combo_i, (p, q) in enumerate(itertools.product(P_GRID, Q_GRID)):
        log(f"--- combinacion {combo_i + 1}/16: p={p}, q={q} ---")
        tc0 = time.time()
        embeddings, timings, _ = train_node2vec(
            args.train_graph, p, q, GRID_WALK_LENGTH, GRID_NUM_WALKS, dim=DIM,
            workers=args.workers, seed=SEED, indptr=indptr, indices=indices, nodes=nodes,
        )
        log(f"  caminatas: {timings['t_walks']:.1f}s ({timings['n_walks']:,} caminatas x "
            f"{timings['walk_length']} pasos) | word2vec: {timings['t_word2vec']:.1f}s | "
            f"embeddings faltantes: {timings['n_missing_embeddings']}")

        for task in TASKS:
            train_pairs, train_y, val_pairs, val_y = task_data[task]
            X_train, kept_train = build_features(embeddings, train_pairs, missing_policy="zero")
            X_val, kept_val = build_features(embeddings, val_pairs, missing_policy="zero")
            _, diag = train_classifier(X_train, train_y, X_val, val_y)
            val_auprc = diag["val_auprc_mejor_c"]
            results.append({
                "p": p, "q": q, "task": task, "val_auprc": val_auprc,
                "best_c": diag["best_c"],
            })
            log(f"    {task}: val_auprc={val_auprc:.4f} (C={diag['best_c']})")

        log(f"  combinacion {combo_i + 1}/16 completa en {time.time() - tc0:.1f}s "
            f"(acumulado {time.time() - t0:.1f}s)")

    df = pd.DataFrame(results)
    df.to_csv(f"{args.out_dir}/grid_results.tsv", sep="\t", index=False)
    log(f"Grid completa -> {args.out_dir}/grid_results.tsv")

    best_per_task = df.loc[df.groupby("task")["val_auprc"].idxmax()]
    best_per_task = best_per_task[["task", "p", "q", "val_auprc", "best_c"]].reset_index(drop=True)
    best_per_task.to_csv(f"{args.out_dir}/best_pq_per_task.tsv", sep="\t", index=False)
    log("Mejor (p,q) por tarea:")
    for _, row in best_per_task.iterrows():
        log(f"  {row['task']}: p={row['p']}, q={row['q']}, val_auprc={row['val_auprc']:.4f}")

    with open(f"{args.out_dir}/grid_manifest.json", "w") as f:
        json.dump({
            "p_grid": P_GRID, "q_grid": Q_GRID, "dim": DIM,
            "grid_walk_length": GRID_WALK_LENGTH, "grid_num_walks": GRID_NUM_WALKS,
            "seed": SEED, "n_combinaciones": len(P_GRID) * len(Q_GRID),
            "tiempo_total_s": time.time() - t0,
            "nota": "Presupuesto reducido usado SOLO para rankear p,q por val_auprc. "
                    "La configuracion final gana con num_walks=10/walk_length=80.",
        }, f, indent=2, ensure_ascii=False)

    log(f"Tiempo total grid: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
