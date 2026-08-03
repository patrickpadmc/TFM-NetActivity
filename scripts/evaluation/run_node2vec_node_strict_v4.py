#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): entrenamiento de Node2Vec a presupuesto
completo (num_walks=10, walk_length=80, igual que run_node2vec_final_v4.py)
para la config (p,q) ganadora de GEN-ass-DIS segun la grid, evaluado en
el split node_strict (unica tarea con este split). Independiente del
job de node2vec_final (no reutiliza ni espera su resultado) porque ese
script no guarda los embeddings en disco -- se reentrena una sola vez
esta config puntual (p=2.0, q=0.25 segun best_pq_per_task.tsv) para
poder evaluar node_strict sin bloquear el resto del pipeline.

Uso:
    python3 run_node2vec_node_strict_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --p 2.0 --q 0.25 \
        --out data/processed/evaluation/graph_v4/node2vec_final/node2vec_results_node_strict.tsv
"""
import argparse
import time

import numpy as np
import pandas as pd

from node2vec_common import train_node2vec
from common_eval import build_features, train_classifier, evaluate

FULL_WALK_LENGTH = 80
FULL_NUM_WALKS = 10
DIM = 128
SEED = 42
TASK = "GEN-ass-DIS"


def log(msg):
    print(f"[n2v_node_strict] {msg}", flush=True)


def load_pairs_labels(task_dir, split, strict=False):
    suffix = "_strict" if strict else ""
    pos = pd.read_csv(f"{task_dir}/{split}{suffix}.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}{suffix}.tsv", sep="\t", dtype=str)
    pairs = list(zip(pos["n1"], pos["n2"])) + list(zip(neg["n1"], neg["n2"]))
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    return pairs, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-graph", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--p", type=float, required=True)
    ap.add_argument("--q", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    t0 = time.time()
    log(f"Entrenando Node2Vec a presupuesto completo (num_walks={FULL_NUM_WALKS}, "
        f"walk_length={FULL_WALK_LENGTH}) con p={args.p}, q={args.q} (config ganadora de {TASK})")
    embeddings, timings, _ = train_node2vec(
        args.train_graph, args.p, args.q, FULL_WALK_LENGTH, FULL_NUM_WALKS, dim=DIM,
        workers=args.workers, seed=SEED,
    )
    log(f"caminatas: {timings['t_walks']:.1f}s ({timings['n_walks']:,} caminatas x "
        f"{timings['walk_length']} pasos) | word2vec: {timings['t_word2vec']:.1f}s | "
        f"total: {time.time() - t0:.1f}s")

    task_dir = f"{args.eval_dir}/{TASK}"
    train_pairs, train_y = load_pairs_labels(task_dir, "train", strict=True)
    val_pairs, val_y = load_pairs_labels(task_dir, "val", strict=True)
    test_pairs, test_y = load_pairs_labels(task_dir, "test", strict=True)

    X_train, _ = build_features(embeddings, train_pairs, missing_policy="zero")
    X_val, _ = build_features(embeddings, val_pairs, missing_policy="zero")
    X_test, _ = build_features(embeddings, test_pairs, missing_policy="zero")

    clf, diag = train_classifier(X_train, train_y, X_val, val_y)
    metrics = evaluate(clf, X_test, test_y)

    result = {
        "task": TASK, "split_type": "node_strict", "method": "node2vec",
        "p": args.p, "q": args.q,
        "auprc": metrics["auprc"], "auroc": metrics["auroc"],
        "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
        "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
    }
    log(f"  {TASK} (node_strict): test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} "
        f"(C={diag['best_c']})")

    pd.DataFrame([result]).to_csv(args.out, sep="\t", index=False)
    log(f"Resultado -> {args.out}")
    log(f"Tiempo total: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
