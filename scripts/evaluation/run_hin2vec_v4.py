#!/usr/bin/env python3
"""
Seccion 5 (TFM graph_v4): entrenamiento de HIN2Vec y evaluacion en las
7 tareas primarias, con el mismo common_eval.py que el resto de metodos.

PRIMER USO (este script, tal cual): PROFILING a escala completa de
graph_v4 -- las 6 relaciones base de los metapaths son intrinsecamente
chicas (10K-315K aristas cada una, muy por debajo de las ~110M del
grafo completo que sí forzaron el uso de un subgrafo para Spectral/GAE/
GAT-AE en la Seccion 4), asi que antes de decidir fracciones de
muestreo (25/50/75/100%) para el subgrafo de benchmark de la Seccion 5,
se mide el costo real de HIN2Vec a escala completa. Si es viable, el
"100%" del diseño multi-escala puede ser el propio graph_v4 (o cerca);
si no, se define el techo empiricamente igual que se hizo con Spectral/
GAE/GAT-AE en la Seccion 4.

Uso:
    python3 run_hin2vec_v4.py \
        --edges data/processed/integrated/graph_v4/graph_edges.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/evaluation/graph_v4/hin2vec \
        --max-epochs 200 --patience 20
"""
import argparse
import json
import os
import resource
import time

import numpy as np
import pandas as pd

from hin2vec_common import train_hin2vec
from common_eval import build_features, train_classifier, evaluate

TASKS = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
         "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]


def log(msg):
    print(f"[run_hin2vec] {msg}", flush=True)


def load_pairs_labels(task_dir, split):
    pos = pd.read_csv(f"{task_dir}/{split}.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}.tsv", sep="\t", dtype=str)
    pairs = list(zip(pos["n1"], pos["n2"])) + list(zip(neg["n1"], neg["n2"]))
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    return pairs, y


def peak_mem_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", required=True)
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--walk-length", type=int, default=6)
    ap.add_argument("--max-epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=20)
    ap.add_argument("--check-every", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    t0 = time.time()
    embeddings, timings = train_hin2vec(
        args.edges, dim=args.dim, walk_length=args.walk_length,
        max_epochs=args.max_epochs, patience=args.patience, check_every=args.check_every,
        seed=args.seed, log_prefix="hin2vec",
    )
    log(f"Entrenamiento HIN2Vec completo. t={time.time() - t0:.1f}s. "
        f"epocas={timings['n_epochs_run']}, mejor val_auprc interno={timings['best_val_auprc_interno']:.4f}, "
        f"memoria pico={peak_mem_gb():.2f}GB")

    np.savez(f"{args.out_dir}/hin2vec_embeddings.npz",
              node_ids=np.array(list(embeddings.keys()), dtype=object),
              embeddings=np.stack(list(embeddings.values())))
    with open(f"{args.out_dir}/hin2vec_training_history.json", "w") as f:
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
            "task": task, "split_type": "edge", "method": "hin2vec",
            "auprc": metrics["auprc"], "auroc": metrics["auroc"],
            "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
            "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
        })
        log(f"  {task}: test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} (C={diag['best_c']})")

    out_df = pd.DataFrame(results)
    out_df.to_csv(f"{args.out_dir}/hin2vec_results.tsv", sep="\t", index=False)
    log(f"Resultados -> {args.out_dir}/hin2vec_results.tsv")
    log(f"Tiempo total: {time.time() - t0:.1f}s. Memoria pico: {peak_mem_gb():.2f}GB")


if __name__ == "__main__":
    main()
