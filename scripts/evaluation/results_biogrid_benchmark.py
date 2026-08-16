from pathlib import Path
import pandas as pd

d = Path("data/processed/analysis/functional_benchmark_v1/biogrid")
files = sorted(d.glob("*__biogrid_results.tsv"))
print(f"\nRESULT_FILES={len(files)}")
for f in files:
    x = pd.read_csv(f, sep="\t")
    print(f"{f.name}\trows={len(x)}\tok={(x.status == 'ok').sum()}")
