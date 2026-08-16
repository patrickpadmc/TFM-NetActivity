
#!/usr/bin/env python3
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score, f1_score, precision_recall_curve,
    precision_score, recall_score, roc_auc_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

CV_SPLITS = ["cv_fold1", "cv_fold2", "cv_fold3"]
C_VALUES = [0.1, 1, 10, 100, 1000]


def model(c):
    return make_pipeline(
        StandardScaler(),
        SVC(C=c, kernel="rbf", class_weight="balanced", cache_size=4096),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-name", required=True)
    ap.add_argument("--source-path", required=True)
    ap.add_argument("--data-path", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{args.source_name}__go_results.tsv"

    go = pd.read_csv(args.data_path, sep="\t", dtype={"ensg": str})
    emb = pd.read_csv(args.source_path, sep="\t")
    ecols = sorted(
        [c for c in emb.columns if c.startswith("e_")],
        key=lambda c: int(c.split("_")[1]),
    )
    emb = emb.set_index("node_id")[ecols]

    rows = []
    for term in sorted(go["go_term"].unique()):
        t0 = time.time()
        x = go[go["go_term"] == term].copy()
        X = emb.loc[x["ensg"]].to_numpy(dtype=np.float32)
        y = x["label"].to_numpy(dtype=int)
        split = x["split"].to_numpy()

        try:
            cv_auc = {}
            for c in C_VALUES:
                fold_scores = []
                for val in CV_SPLITS:
                    train_mask = np.isin(split, [s for s in CV_SPLITS if s != val])
                    val_mask = split == val
                    clf = model(c)
                    clf.fit(X[train_mask], y[train_mask])
                    fold_scores.append(
                        roc_auc_score(y[val_mask], clf.decision_function(X[val_mask]))
                    )
                cv_auc[c] = float(np.mean(fold_scores))

            best_c = max(cv_auc, key=cv_auc.get)

            oof_y, oof_scores = [], []
            for val in CV_SPLITS:
                train_mask = np.isin(split, [s for s in CV_SPLITS if s != val])
                val_mask = split == val
                clf = model(best_c)
                clf.fit(X[train_mask], y[train_mask])
                oof_y.extend(y[val_mask])
                oof_scores.extend(clf.decision_function(X[val_mask]))

            precision, recall, thresholds = precision_recall_curve(oof_y, oof_scores)
            f1_values = 2 * precision[:-1] * recall[:-1] / (
                precision[:-1] + recall[:-1] + 1e-12
            )
            threshold = float(thresholds[np.argmax(f1_values)])
            val_f1 = float(np.max(f1_values))

            train_mask = np.isin(split, CV_SPLITS)
            test_mask = split == "holdout"
            clf = model(best_c)
            clf.fit(X[train_mask], y[train_mask])
            scores = clf.decision_function(X[test_mask])
            pred = scores >= threshold

            row = {
                "source": args.source_name,
                "go_term": term,
                "dimension": len(ecols),
                "best_c": best_c,
                "cv_mean_auroc": cv_auc[best_c],
                "validation_f1_threshold": threshold,
                "validation_f1": val_f1,
                "n_train": int(train_mask.sum()),
                "n_test": int(test_mask.sum()),
                "n_test_positive": int(y[test_mask].sum()),
                "test_auprc": average_precision_score(y[test_mask], scores),
                "test_auroc": roc_auc_score(y[test_mask], scores),
                "test_f1": f1_score(y[test_mask], pred, zero_division=0),
                "test_precision": precision_score(y[test_mask], pred, zero_division=0),
                "test_recall": recall_score(y[test_mask], pred, zero_division=0),
                "elapsed_seconds": round(time.time() - t0, 2),
                "status": "ok",
            }
        except Exception as exc:
            row = {
                "source": args.source_name,
                "go_term": term,
                "dimension": len(ecols),
                "status": "error",
                "error": repr(exc),
                "elapsed_seconds": round(time.time() - t0, 2),
            }

        rows.append(row)
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
        print(f"[{args.source_name}] {term}: {row['status']}", flush=True)

    print(f"Results: {out_path}", flush=True)


if __name__ == "__main__":
    main()
