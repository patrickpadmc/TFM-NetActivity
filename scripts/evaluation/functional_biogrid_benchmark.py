from pathlib import Path
import pandas as pd

out = Path("data/processed/analysis/functional_benchmark_v1")
out.mkdir(parents=True, exist_ok=True)

# GO: media y mediana entre los 56 términos
go_dir = out / "go"
go = pd.concat(
    [pd.read_csv(f, sep="\t") for f in sorted(go_dir.glob("*__go_results.tsv"))],
    ignore_index=True,
)
go = go[go["status"] == "ok"]

go_summary = (
    go.groupby("source")[["test_auprc", "test_auroc", "test_f1"]]
    .agg(["mean", "median"])
    .reset_index()
)
go_summary.columns = [
    "source", "go_auprc_mean", "go_auprc_median",
    "go_auroc_mean", "go_auroc_median",
    "go_f1_mean", "go_f1_median",
]

# BioGRID: media entre tres folds para cada tarea
bg_dir = out / "biogrid"
bg = pd.concat(
    [pd.read_csv(f, sep="\t") for f in sorted(bg_dir.glob("*__biogrid_results.tsv"))],
    ignore_index=True,
)
bg = bg[bg["status"] == "ok"]

bg_summary = (
    bg.groupby(["source", "task"])[["test_auprc", "test_auroc", "test_pr_at_10"]]
    .mean()
    .reset_index()
    .pivot(index="source", columns="task")
)
bg_summary.columns = [
    f"{task}_{metric}"
    for metric, task in bg_summary.columns
]
bg_summary = bg_summary.reset_index()

master = go_summary.merge(bg_summary, on="source", how="outer")
master = master.sort_values("go_auprc_mean", ascending=False)
master.to_csv(out / "functional_benchmark_master_summary.tsv", sep="\t", index=False)

go.to_csv(out / "go_all_term_results.tsv", sep="\t", index=False)
bg.to_csv(out / "biogrid_all_fold_results.tsv", sep="\t", index=False)

print(master.round(4).to_string(index=False))
print(f"\nSaved: {out / 'functional_benchmark_master_summary.tsv'}")
