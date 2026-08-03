#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): entrenamiento FINAL de Node2Vec (presupuesto
completo: num_walks=10, walk_length=80, decision confirmada) para cada
(p,q) ganador por tarea (segun node2vec_grid/best_pq_per_task.tsv),
seguido de evaluacion en TEST una unica vez por tarea, tal como exige
el protocolo ("Ejecuta el test una unica vez con la mejor
configuracion").
Varias tareas pueden compartir el mismo (p,q) ganador -- en ese caso
se entrena una sola vez ese embedding y se reutiliza para todas.

Seccion 9: Node2Vec no produce UN embedding global -- entrena 5
configuraciones (p,q) distintas, una por grupo de tareas. Para el
catalogo unificado de embeddings se persiste como representativo el
embedding de la configuracion (p,q) que incluye GEN-ass-DIS (tarea
insignia usada en Secciones 6-8), documentado explicitamente como
simplificacion metodologica -- Node2Vec es inherentemente especifico
por tarea, esto no es "el" embedding de Node2Vec en sentido estricto.

Uso:
    python3 run_node2vec_final_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --best-pq data/processed/evaluation/graph_v4/node2vec_grid/best_pq_per_task.tsv \
        --out-dir data/processed/evaluation/graph_v4/node2vec_final
"""
import argparse
import json
import os
import time
import numpy as np
import pandas as pd
from node2vec_common import build_csr, train_node2vec
from common_eval import build_features, train_classifier, evaluate

FULL_WALK_LENGTH = 80
FULL_NUM_WALKS = 10
DIM = 128
SEED = 42
REPRESENTATIVE_TASK = "GEN-ass-DIS"


def log(msg):
    print(f"[n2v_final] {msg}", flush=True)


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
    ap.add_argument("--best-pq", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    SEED = args.seed

    os.makedirs(args.out_dir, exist_ok=True)

    best_pq = pd.read_csv(args.best_pq, sep="\t")
    log("Mejores (p,q) por tarea (de la grid):")
    for _, row in best_pq.iterrows():
        log(f"  {row['task']}: p={row['p']}, q={row['q']}, val_auprc(grid)={row['val_auprc']:.4f}")

    distinct_pq = best_pq[["p", "q"]].drop_duplicates().to_records(index=False).tolist()
    log(f"Configuraciones (p,q) distintas a entrenar a presupuesto completo: {len(distinct_pq)}")

    t0 = time.time()
    indptr, indices, nodes, _ = build_csr(args.train_graph)
    log(f"{len(nodes):,} nodos. CSR listo en {time.time() - t0:.1f}s")

    results = []
    representative_saved = False
    for p, q in distinct_pq:
        tasks_for_this_pq = best_pq.loc[(best_pq["p"] == p) & (best_pq["q"] == q), "task"].tolist()
        log(f"--- (p,q)=({p},{q}) -> tareas: {tasks_for_this_pq} ---")
        tc0 = time.time()
        embeddings, timings, _ = train_node2vec(
            args.train_graph, p, q, FULL_WALK_LENGTH, FULL_NUM_WALKS, dim=DIM,
            workers=args.workers, seed=SEED, indptr=indptr, indices=indices, nodes=nodes,
        )
        log(f"  caminatas: {timings['t_walks']:.1f}s ({timings['n_walks']:,} caminatas x "
            f"{timings['walk_length']} pasos) | word2vec: {timings['t_word2vec']:.1f}s")

        if REPRESENTATIVE_TASK in tasks_for_this_pq:
            np.savez(f"{args.out_dir}/node2vec_embeddings.npz",
                     node_ids=np.array(list(embeddings.keys()), dtype=object),
                     embeddings=np.stack(list(embeddings.values())))
            with open(f"{args.out_dir}/node2vec_embeddings_representative_config.json", "w") as f:
                json.dump({"p": p, "q": q, "seed": SEED,
                           "tarea_representativa": REPRESENTATIVE_TASK,
                           "tareas_que_comparten_esta_config": tasks_for_this_pq,
                           "nota": "Node2Vec entrena 5 configuraciones (p,q) distintas; esta es la "
                                   "unica persistida como representativa del catalogo unificado "
                                   "de la Seccion 9, por ser la config ganadora de GEN-ass-DIS."},
                          f, indent=2, ensure_ascii=False)
            representative_saved = True
            log(f"  >> Config representativa para el catalogo (incluye {REPRESENTATIVE_TASK}) -> node2vec_embeddings.npz")

        for task in tasks_for_this_pq:
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
                "task": task, "split_type": "edge", "method": "node2vec",
                "p": p, "q": q, "auprc": metrics["auprc"], "auroc": metrics["auroc"],
                "n_test": metrics["n_test"], "n_pos_test": metrics["n_pos_test"],
                "best_c": diag["best_c"], "val_auprc_mejor_c": diag["val_auprc_mejor_c"],
            })
            log(f"    {task}: test_auprc={metrics['auprc']:.4f} test_auroc={metrics['auroc']:.4f} (C={diag['best_c']})")
        log(f"  (p,q)=({p},{q}) completo en {time.time() - tc0:.1f}s (acumulado {time.time() - t0:.1f}s)")

    if not representative_saved:
        log(f"AVISO: {REPRESENTATIVE_TASK} no aparecio en ninguna config -- no se guardo ningun .npz representativo")

    out_df = pd.DataFrame(results)
    out_df.to_csv(f"{args.out_dir}/node2vec_final_results.tsv", sep="\t", index=False)
    log(f"Resultados -> {args.out_dir}/node2vec_final_results.tsv")

    with open(f"{args.out_dir}/node2vec_final_manifest.json", "w") as f:
        json.dump({
            "walk_length": FULL_WALK_LENGTH, "num_walks": FULL_NUM_WALKS, "dim": DIM, "seed": SEED,
            "n_configuraciones_pq_distintas": len(distinct_pq),
            "tiempo_total_s": time.time() - t0,
        }, f, indent=2, ensure_ascii=False)
    log(f"Tiempo total: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
