#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, make_scorer, roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC

C_VALUES = [0.1, 1, 10, 100, 1000]


def precision_at_k(y_true, scores, k=10):
    order = np.argsort(scores)[::-1]
    return float(np.mean(np.asarray(y_true)[order[:k]]))


def pair_features(df, embeddings):
    return (
        embeddings.loc[df["ensg_a"]].to_numpy(dtype=np.float32)
        + embeddings.loc[df["ensg_b"]].to_numpy(dtype=np.float32)
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-name", required=True)
    ap.add_argument("--source-path", required=True)
    ap.add_argument("--splits-path", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--task", choices=["synthetic_lethality", "negative_genetic"])
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    emb = pd.read_csv(args.source_path, sep="\t")
    cols = sorted(
        [c for c in emb.columns if c.startswith("e_")],
        key=lambda c: int(c.split("_")[1]),
    )
    embeddings = emb.set_index("node_id")[cols]

    pairs = pd.read_csv(args.splits_path, sep="\t", dtype={
        "ensg_a": str, "ensg_b": str, "pair_index_reference": int,
    })
    pairs["label"] = pairs["label"].astype(int)
    tasks = [args.task] if args.task else sorted(pairs["task"].unique())

    all_rows = []
    for task in tasks:
        for outer_fold in sorted(pairs["outer_fold"].unique()):
            outer = pairs[
                (pairs["task"] == task)
                & (pairs["outer_fold"] == outer_fold)
            ]
            train = outer[outer["partition"] == "outer_train"].copy()
            test = outer[outer["partition"] == "outer_test"].copy()

            X_train, y_train = pair_features(train, embeddings), train["label"].to_numpy()
            X_test, y_test = pair_features(test, embeddings), test["label"].to_numpy()

            position = {
                idx: pos for pos, idx in enumerate(train["pair_index_reference"])
            }
            inner_cv = []
            for inner_fold in [0, 1, 2]:
                itr = outer[
                    (outer["partition"] == "inner_train")
                    & (outer["inner_fold"] == inner_fold)
                ]["pair_index_reference"]
                iva = outer[
                    (outer["partition"] == "inner_validation")
                    & (outer["inner_fold"] == inner_fold)
                ]["pair_index_reference"]
                inner_cv.append((
                    [position[i] for i in itr if i in position],
                    [position[i] for i in iva if i in position],
                ))

            scorer = make_scorer(precision_at_k, response_method="decision_function", k=10)
            grid = GridSearchCV(
                SVC(kernel="rbf", class_weight="balanced", cache_size=4096),
                param_grid={"C": C_VALUES},
                scoring={"auroc": "roc_auc", "auprc": "average_precision", "pr_at_10": scorer},
                refit="auroc",
                cv=inner_cv,
                n_jobs=3,
                return_train_score=False,
            )
            grid.fit(X_train, y_train)

            scores = grid.decision_function(X_test)
            best = grid.best_index_
            result = {
                "source": args.source_name,
                "task": task,
                "outer_fold": outer_fold,
                "dimension": len(cols),
                "best_c": grid.best_params_["C"],
                "inner_auroc": grid.best_score_,
                "inner_auprc": grid.cv_results_["mean_test_auprc"][best],
                "inner_pr_at_10": grid.cv_results_["mean_test_pr_at_10"][best],
                "n_train": len(y_train),
                "n_test": len(y_test),
                "n_test_positive": int(y_test.sum()),
                "test_auroc": roc_auc_score(y_test, scores),
                "test_auprc": average_precision_score(y_test, scores),
                "test_pr_at_10": precision_at_k(y_test, scores),
                "status": "ok",
            }
            all_rows.append(result)
            pd.DataFrame(all_rows).to_csv(
                out_dir / f"{args.source_name}__biogrid_results.tsv",
                sep="\t", index=False,
            )
            print(f"[{args.source_name}] {task} fold={outer_fold}: ok", flush=True)


if __name__ == "__main__":
    main()
