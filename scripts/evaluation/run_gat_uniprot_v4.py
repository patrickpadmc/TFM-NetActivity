#!/usr/bin/env python3
"""
Seccion 10 (TFM graph_v4): entrenamiento de GAT-AE sobre el subgrafo
comun con inicializacion de nodos GEN via UniProt (condiciones
uniprot_fixed / uniprot_finetuned), evaluado en las 7 tareas primarias
con el mismo common_eval.py que el resto de metodos. Analogo a
run_gat_v4.py (condicion 'structural'), pero pasando
--uniprot-tsv/--uniprot-mode a train_gat (ver gat_common.py).

Persiste, ademas de lo que ya guarda run_gat_v4.py:
  - gat_uniprot_{mode}_initial_embeddings.npz: vectores UniProt
    proyectados ANTES de entrenar (solo para los nodos GEN cubiertos),
    exigido por el protocolo de la Seccion 10 junto a los finales.
  - gat_training_history.json: igual que en run_gat_v4.py, pero SIN
    los propios vectores iniciales (no son serializables a JSON de
    forma directa -- van en el .npz de arriba).

Uso:
    python3 run_gat_uniprot_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/evaluation/graph_v4/gat_uniprot_fixed \
        --uniprot-tsv data/processed/analysis/external_embeddings_v4/uniprot_embeddings_subgraph.tsv \
        --uniprot-mode fixed \
        --weight-decay 0.0
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
    print(f"[run_gat_uniprot] {msg}", flush=True)


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
    ap.add_argument("--uniprot-tsv", required=True)
    ap.add_argument("--uniprot-mode", required=True, choices=["fixed", "finetuned"])
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--max-epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=10)
    ap.add_argument("--check-every", type=int, default=5)
    ap.add_argument("--weight-decay", type=float, default=5e-4)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    t0 = time.time()
    embeddings, timings = train_gat(
        args.train_graph, dim=args.dim,
        max_epochs=args.max_epochs, patience=args.patience, check_every=args.check_every,
        weight_decay=args.weight_decay, seed=args.seed,
        log_prefix=f"gat_uniprot_{args.uniprot_mode}",
        uniprot_tsv=args.uniprot_tsv, uniprot_mode=args.uniprot_mode,
    )
    log(f"Entrenamiento GAT-AE (uniprot_{args.uniprot_mode}) completo. t={time.time() - t0:.1f}s. "
        f"epocas={timings['n_epochs_run']}, mejor val_auprc interno={timings['best_val_auprc_interno']:.4f}, "
        f"nodos GEN cubiertos por UniProt={timings['n_uniprot_covered']:,}")

    np.savez(f"{args.out_dir}/gat_embeddings.npz",
              node_ids=np.array(list(embeddings.keys()), dtype=object),
              embeddings=np.stack(list(embeddings.values())))

    # Vectores iniciales (UniProt proyectado, ANTES de entrenar) -- exigido
    # junto a los finales por el protocolo de la Seccion 10. Se guardan
    # aparte del .npz de embeddings finales porque cubren solo los nodos
    # GEN con cobertura UniProt, no todos los nodos del grafo.
    initial_emb = timings.pop("initial_uniprot_embeddings")
    if initial_emb:
        np.savez(f"{args.out_dir}/gat_uniprot_{args.uniprot_mode}_initial_embeddings.npz",
                  node_ids=np.array(list(initial_emb.keys()), dtype=object),
                  embeddings=np.stack(list(initial_emb.values())))
        log(f"Vectores iniciales UniProt proyectados guardados "
            f"({len(initial_emb):,} nodos) -> gat_uniprot_{args.uniprot_mode}_initial_embeddings.npz")

    # timings ya no contiene arrays de numpy tras el pop() de arriba -> JSON-safe
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
            "task": task, "split_type": "edge", "method": f"gatae_uniprot_{args.uniprot_mode}",
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
