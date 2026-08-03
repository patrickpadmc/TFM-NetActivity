#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): entrenamiento de GAT-AE sobre el subgrafo
comun (ver docstring de gat_common.py para la justificacion de memoria)
y evaluacion en las 7 tareas primarias, con el mismo common_eval.py
que el resto de metodos.

Uso:
    python3 run_gat_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/evaluation/graph_v4/gat
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd

from gat_common import train_gat
from common_eval import build_features, train_classifier, evaluate

TASKS = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
         "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]


def log(msg):
    print(f"[run_gat] {msg}", flush=True)


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
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--max-epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=10)
    ap.add_argument("--check-every", type=int, default=5)
    ap.add_argument("--weight-decay", type=float, default=5e-4)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    t0 = time.time()
    embeddings, timings = train_gat(
        args.train_graph, dim=args.dim,
        max_epochs=args.max_epochs, patience=args.patience, check_every=args.check_every,
        weight_decay=args.weight_decay, log_prefix="gat",
    )
    log(f"Entrenamiento GAT-AE completo. t={time.time() - t0:.1f}s. "
        f"epocas={timings['n_epochs_run']}, mejor val_auprc interno={timings['best_val_auprc_interno']:.4f}")

    np.savez(f"{args.out_dir}/gat_embeddings.npz",
              node_ids=np.array(list(embeddings.keys()), dtype=object),
              embeddings=np.stack(list(embeddings.values())))
    with open(f"{args.out_dir}/gat_training_history.json", "w") as f:
        json.dump(timings, f, indent=2, ensure_ascii=False)

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
            "task": task, "split_type": "edge", "method": "gatae",
            "auprc": metrics["auprc"], "auroc": metrics["auroc"],
            "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
            "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
        })
        log(f"  {task}: test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} (C={diag['best_c']})")

    out_df = pd.DataFrame(results)
    out_df.to_csv(f"{args.out_dir}/gat_results.tsv", sep="\t", index=False)
    log(f"Resultados -> {args.out_dir}/gat_results.tsv")
    log(f"Tiempo total (entrenamiento + evaluacion 7 tareas): {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
