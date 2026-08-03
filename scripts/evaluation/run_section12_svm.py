#!/usr/bin/env python3
"""
Seccion 12 (TFM graph_v4): predice enlaces con SVM (lineal + RBF
aproximado por Nystroem) para UNA fuente de embedding, iterando sobre
todas las relaciones aplicables a esa fuente (7 relaciones para fuentes
de grafo completo; solo GEN-ppi-GEN para UniProt/GenePT crudos, que
solo cubren nodos GEN).

Reutiliza EXACTAMENTE los mismos splits train/val/test y pares
negativos ya usados por run_gae_v4.py / run_gat_v4.py (Secciones 9-10),
vía data/processed/evaluation/graph_v4/<relacion>/{split}.tsv +
negatives_{split}.tsv -- variante NO estricta (mismo convenio que el
resto del catalogo, para comparabilidad de resultados entre secciones).

Uso (fuente de grafo completo, ej. GAE estructural):
    python3 run_section12_svm.py \
        --source-name gae \
        --source-path data/processed/analysis/embedding_catalog_v4/gae_structural_seed42.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/analysis/seccion12_svm_link_prediction

Uso (fuente cruda, solo GEN-GEN):
    python3 run_section12_svm.py \
        --source-name uniprot_raw \
        --source-path data/processed/analysis/external_embeddings_v4/uniprot_embeddings_subgraph.tsv \
        --eval-dir data/processed/evaluation/graph_v4 \
        --out-dir data/processed/analysis/seccion12_svm_link_prediction \
        --raw-external
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from svm_link_pred_common import (
    build_pair_features, class_weight_setting, select_best_model,
    select_threshold, evaluate_test, pr_roc_curves, ranked_evaluation_svm, log,
)

TASKS_ALL = ["GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
             "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN"]
RELATION_TYPES = {
    "GEN-ass-DIS": ("GEN", "DIS"), "GEN-ass-TIS": ("GEN", "TIS"),
    "GEN-ass-PWY": ("GEN", "PWY"), "CPD-trt-DIS": ("CPD", "DIS"),
    "CPD-int-GEN": ("CPD", "GEN"), "GEN-ppi-GEN": ("GEN", "GEN"),
    "CLL-mut-GEN": ("CLL", "GEN"),
}
RAW_EXTERNAL_TASKS = ["GEN-ppi-GEN"]
K_RANK = 10
N_NEG_PER_QUERY = 99
SEED = 42


def load_embeddings_full(path):
    """Carga TODOS los tipos de nodo (a diferencia de Seccion 11, aqui
    hacen falta -- las relaciones conectan tipos distintos)."""
    df = pd.read_csv(path, sep="\t")
    e_cols = sorted((c for c in df.columns if c.startswith("e_")), key=lambda c: int(c.split("_")[1]))
    ids = df["node_id"].to_numpy()
    types = df["node_type"].to_numpy()
    mat = df[e_cols].to_numpy(dtype=np.float32)
    embeddings = dict(zip(ids, mat))
    node_types = dict(zip(ids, types))
    return embeddings, node_types, len(e_cols)


def load_pairs_labels(task_dir, split):
    pos = pd.read_csv(f"{task_dir}/{split}.tsv", sep="\t", dtype=str)
    neg = pd.read_csv(f"{task_dir}/negatives_{split}.tsv", sep="\t", dtype=str)
    pos_pairs = list(zip(pos["n1"], pos["n2"]))
    neg_pairs = list(zip(neg["n1"], neg["n2"]))
    pairs = pos_pairs + neg_pairs
    y = np.concatenate([np.ones(len(pos_pairs)), np.zeros(len(neg_pairs))])
    return pairs, y, pos_pairs


def run_task(source_name, task, embeddings, node_types, eval_dir, out_dir, seed):
    t0 = time.time()
    task_dir = f"{eval_dir}/{task}"
    log(f"[{source_name}/{task}] cargando pares...")
    train_pairs, y_train_full, train_pos = load_pairs_labels(task_dir, "train")
    val_pairs, y_val_full, val_pos = load_pairs_labels(task_dir, "val")
    test_pairs, y_test_full, test_pos = load_pairs_labels(task_dir, "test")

    X_train, kept_tr = build_pair_features(embeddings, train_pairs)
    y_train = y_train_full[kept_tr]
    X_val, kept_va = build_pair_features(embeddings, val_pairs)
    y_val = y_val_full[kept_va]
    X_test, kept_te = build_pair_features(embeddings, test_pairs)
    y_test = y_test_full[kept_te]

    n_dropped = dict(
        train=int(len(train_pairs) - len(y_train)),
        val=int(len(val_pairs) - len(y_val)),
        test=int(len(test_pairs) - len(y_test)),
    )
    log(f"[{source_name}/{task}] train={len(y_train)} (desc.{n_dropped['train']}) "
        f"val={len(y_val)} (desc.{n_dropped['val']}) test={len(y_test)} (desc.{n_dropped['test']})")

    if len(y_train) < 20 or len(np.unique(y_train)) < 2 or len(y_val) < 5 or len(np.unique(y_val)) < 2:
        log(f"[{source_name}/{task}] AVISO: datos insuficientes tras el filtrado, se omite esta combinacion.")
        return None

    scaler = StandardScaler().fit(X_train)
    X_train_s = scaler.transform(X_train).astype(np.float32)
    X_val_s = scaler.transform(X_val).astype(np.float32)
    X_test_s = scaler.transform(X_test).astype(np.float32)

    cw, frac_pos_train = class_weight_setting(y_train)
    log(f"[{source_name}/{task}] frac_pos_train={frac_pos_train:.3f} class_weight={cw}")

    model, diag = select_best_model(X_train_s, y_train, X_val_s, y_val, cw=cw, seed=seed)
    log(f"[{source_name}/{task}] modelo elegido={diag['elegido']} "
        f"(lineal val_auprc={diag['linear']['val_auprc']:.4f}, "
        f"rbf_approx val_auprc={diag['rbf_approx']['val_auprc']:.4f})")

    val_scores = model.decision_function(X_val_s)
    threshold, thr_diag = select_threshold(y_val, val_scores)

    test_scores = model.decision_function(X_test_s)
    test_metrics = evaluate_test(y_test, test_scores, threshold)
    log(f"[{source_name}/{task}] TEST auprc={test_metrics['auprc']:.4f} auroc={test_metrics['auroc']:.4f} "
        f"f1={test_metrics['f1']:.4f}")

    curves = pr_roc_curves(y_test, test_scores)
    np.savez(f"{out_dir}/{source_name}__{task}__curves.npz", **curves)

    # --- ranking (Hits@10 / MRR) ---
    t1, t2 = RELATION_TYPES[task]
    node_pool_t2 = np.array([nid for nid, t in node_types.items() if t == t2 and nid in embeddings])
    forbidden = set(train_pos) | set(val_pos) | set(test_pos)
    test_positive_pairs = [p for p, y in zip(test_pairs, y_test_full) if y == 1]
    rng = np.random.default_rng(seed)
    if len(node_pool_t2) >= N_NEG_PER_QUERY + 1 and test_positive_pairs:
        ranked = ranked_evaluation_svm(model, embeddings, test_positive_pairs, node_pool_t2, rng,
                                        k=K_RANK, n_negatives_per_query=N_NEG_PER_QUERY,
                                        forbidden_pairs=forbidden)
    else:
        ranked = {"hits_at_k": None, "mrr": None, "n_queries": 0, "k": K_RANK,
                  "nota": "pool de candidatos t2 insuficiente para ranking"}
    log(f"[{source_name}/{task}] Hits@{K_RANK}={ranked['hits_at_k']} MRR={ranked['mrr']} "
        f"(n_queries={ranked['n_queries']})")

    elapsed = time.time() - t0
    result = dict(
        source=source_name, relation=task,
        n_train=len(y_train), n_val=len(y_val), n_test=len(y_test),
        n_dropped_train=n_dropped["train"], n_dropped_val=n_dropped["val"], n_dropped_test=n_dropped["test"],
        frac_pos_train=frac_pos_train, class_weight=str(cw),
        modelo_elegido=diag["elegido"],
        linear_best_c=diag["linear"]["best_c"], linear_val_auprc=diag["linear"]["val_auprc"],
        rbf_best_c=diag["rbf_approx"]["best_c"], rbf_best_gamma=diag["rbf_approx"]["best_gamma"],
        rbf_n_components=diag["rbf_approx"]["n_components"], rbf_val_auprc=diag["rbf_approx"]["val_auprc"],
        threshold=threshold, f1_val_en_umbral=thr_diag["f1_val"],
        test_auprc=test_metrics["auprc"], test_auroc=test_metrics["auroc"],
        test_f1=test_metrics["f1"], test_precision=test_metrics["precision"], test_recall=test_metrics["recall"],
        hits_at_10=ranked["hits_at_k"], mrr=ranked["mrr"], n_ranked_queries=ranked["n_queries"],
        elapsed_seconds=round(elapsed, 1),
    )
    with open(f"{out_dir}/{source_name}__{task}__result.json", "w") as f:
        json.dump(result, f, indent=2, default=str)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-name", required=True)
    ap.add_argument("--source-path", required=True)
    ap.add_argument("--eval-dir", default="data/processed/evaluation/graph_v4")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--raw-external", action="store_true",
                     help="Solo evalua GEN-ppi-GEN (embeddings crudos UniProt/GenePT, sin cobertura completa de nodos)")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    log(f"Cargando embeddings de {args.source_name} desde {args.source_path}...")
    embeddings, node_types, dim = load_embeddings_full(args.source_path)
    log(f"  {len(embeddings):,} nodos, dim={dim}")

    tasks = RAW_EXTERNAL_TASKS if args.raw_external else TASKS_ALL
    log(f"Relaciones a evaluar: {tasks}")

    rows = []
    for task in tasks:
        res = run_task(args.source_name, task, embeddings, node_types, args.eval_dir, args.out_dir, args.seed)
        if res is not None:
            rows.append(res)

    if rows:
        out_df = pd.DataFrame(rows)
        out_path = f"{args.out_dir}/{args.source_name}__svm_results.tsv"
        out_df.to_csv(out_path, sep="\t", index=False)
        log(f"Tabla de resultados -> {out_path}")
    log("Listo.")


if __name__ == "__main__":
    main()
