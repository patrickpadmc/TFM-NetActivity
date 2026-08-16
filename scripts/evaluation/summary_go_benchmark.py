from pathlib import Path
import pandas as pd

d = Path("data/processed/analysis/functional_benchmark_v1/go")
files = sorted(d.glob("*__go_results.tsv"))
df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)

print(f"RESULT_FILES={len(files)}")
print(f"TOTAL_ROWS={len(df)}")
print(f"OK_ROWS={(df.status == 'ok').sum()}")
print(f"EXPECTED_ROWS={11 * 56}")
print(f"UNIQUE_EMBEDDINGS={df.source.nunique()}")
print(f"UNIQUE_GO_TERMS={df.go_term.nunique()}")

summary = (df[df.status == "ok"]
           .groupby("source")[["test_auprc", "test_auroc", "test_f1"]]
           .agg(["mean", "median"])
           .round(4))
print("\nSummary by embedding:")
print(summary.to_string())
