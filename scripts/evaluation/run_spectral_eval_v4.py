#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): evaluacion de los embeddings de Spectral
Embedding (ya calculados sobre el subgrafo comun, job 308480) en las 7
tareas primarias, con el mismo common_eval.py que el resto de metodos.

Uso:
    python3 run_spectral_eval_v4.py \
        --embeddings data/processed/evaluation/graph_v4/embeddings/spectral_v4_subgraph.npz \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out data/processed/evaluation/graph_v4/spectral_results.tsv
"""
import argparse

import numpy as np
import pandas as pd

from common_eval import build_features, train_classifier, evaluate

TASKS = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
         "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]


def log(msg):
    print(f"[spectral_eval] {msg}", flush=True)


def load_pairs_labels(task_dir, split):
    pos = pd.read_csv(f"{task_dir}/{split}.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}.tsv", sep="\t", dtype=str)
    pairs = list(zip(pos["n1"], pos["n2"])) + list(zip(neg["n1"], neg["n2"]))
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    return pairs, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    log(f"Cargando embeddings -> {args.embeddings}")
    data = np.load(args.embeddings, allow_pickle=True)
    node_ids = data["node_ids"]
    vecs = data["embeddings"]
    embeddings = {node_ids[i]: vecs[i] for i in range(len(node_ids))}
    log(f"{len(embeddings):,} nodos con embedding (dim={vecs.shape[1]})")

    results = []
    for task in TASKS:
        task_dir = f"{args.eval_dir}/{task}"
        train_pairs, train_y = load_pairs_labels(task_dir, "train")
        val_pairs, val_y = load_pairs_labels(task_dir, "val")
        test_pairs, test_y = load_pairs_labels(task_dir, "test")

        X_train, _ = build_features(embeddings, train_pairs, missing_policy="zero")
        X_val, _ = build_features(embeddings, val_pairs, missing_policy="zero")
        X_test, _ = build_features(embeddings, test_pairs, missing_policy="zero")

        clf, diag = train_classifier(X_train, train_y, X_val, val_y)
        metrics = evaluate(clf, X_test, test_y)

        results.append({
            "task": task, "split_type": "edge", "method": "spectral",
            "auprc": metrics["auprc"], "auroc": metrics["auroc"],
            "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
            "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
        })
        log(f"  {task}: test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} (C={diag['best_c']})")

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, sep="\t", index=False)
    log(f"Resultados -> {args.out}")


if __name__ == "__main__":
    main()
