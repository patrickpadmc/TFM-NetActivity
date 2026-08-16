from pathlib import Path
import pandas as pd

d = Path("data/processed/analysis/functional_benchmark_v1/biogrid")
df = pd.concat(
    [pd.read_csv(f, sep="\t") for f in sorted(d.glob("*__biogrid_results.tsv"))],
    ignore_index=True,
)

print(f"TOTAL_ROWS={len(df)}")
print(f"OK_ROWS={(df.status == 'ok').sum()}")

summary = (
    df[df.status == "ok"]
    .groupby(["task", "source"])[["test_auprc", "test_auroc", "test_pr_at_10"]]
    .mean()
    .round(4)
)
print("\nMean metrics by task and embedding:")
print(summary.to_string())

overall = (
    df[df.status == "ok"]
    .groupby("source")[["test_auprc", "test_auroc", "test_pr_at_10"]]
    .mean()
    .sort_values("test_auprc", ascending=False)
    .round(4)
)
print("\nOverall mean across SL + NG:")
print(overall.to_string())
