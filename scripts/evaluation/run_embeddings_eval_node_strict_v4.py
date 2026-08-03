#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): evaluacion generica en el split node_strict de
GEN-ass-DIS (unica tarea con este split) para cualquier metodo cuyos
embeddings ya esten guardados en un .npz (node_ids, embeddings) --
spectral, GAE, GAT-AE. No reentrena nada, solo evalua sobre los pares
*_strict.tsv reutilizando el mismo embedding ya usado en el split "edge".

Uso:
    python3 run_embeddings_eval_node_strict_v4.py \
        --embeddings data/processed/evaluation/graph_v4/gae/gae_embeddings.npz \
        --method gae \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out data/processed/evaluation/graph_v4/gae/gae_results_node_strict.tsv
"""
import argparse

import numpy as np
import pandas as pd

from common_eval import build_features, train_classifier, evaluate

TASK = "GEN-ass-DIS"


def load_pairs_labels(task_dir, split):
    pos = pd.read_csv(f"{task_dir}/{split}_strict.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}_strict.tsv", sep="\t", dtype=str)
    pairs = list(zip(pos["n1"], pos["n2"])) + list(zip(neg["n1"], neg["n2"]))
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    return pairs, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    def log(msg):
        print(f"[eval_node_strict:{args.method}] {msg}", flush=True)

    log(f"Cargando embeddings -> {args.embeddings}")
    data = np.load(args.embeddings, allow_pickle=True)
    node_ids = data["node_ids"]
    vecs = data["embeddings"]
    embeddings = {node_ids[i]: vecs[i] for i in range(len(node_ids))}
    log(f"{len(embeddings):,} nodos con embedding (dim={vecs.shape[1]})")

    task_dir = f"{args.eval_dir}/{TASK}"
    train_pairs, train_y = load_pairs_labels(task_dir, "train")
    val_pairs, val_y = load_pairs_labels(task_dir, "val")
    test_pairs, test_y = load_pairs_labels(task_dir, "test")

    X_train, _ = build_features(embeddings, train_pairs, missing_policy="zero")
    X_val, _ = build_features(embeddings, val_pairs, missing_policy="zero")
    X_test, _ = build_features(embeddings, test_pairs, missing_policy="zero")

    clf, diag = train_classifier(X_train, train_y, X_val, val_y)
    metrics = evaluate(clf, X_test, test_y)

    result = {
        "task": TASK, "split_type": "node_strict", "method": args.method,
        "auprc": metrics["auprc"], "auroc": metrics["auroc"],
        "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
        "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
    }
    log(f"  {TASK} (node_strict): test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} "
        f"(C={diag['best_c']})")

    pd.DataFrame([result]).to_csv(args.out, sep="\t", index=False)
    log(f"Resultado -> {args.out}")


if __name__ == "__main__":
    main()
