#!/usr/bin/env python3
"""
Seccion 6 (TFM graph_v4): construye una variante de ablacion del
subgrafo comun (100%) retirando POR COMPLETO una relacion (familia)
especifica, para medir su contribucion al desempeno de HIN2Vec en
GEN-ass-DIS. No se re-muestrea nada mas: se parte del subgrafo comun
de la Seccion 4/5 y se filtran solo las aristas de la relacion indicada.

Uso:
    python3 build_ablation_subgraph_v4.py \
        --base-subgraph data/processed/evaluation/graph_v4/benchmark_multiscale/graph_benchmark_100pct_v4.tsv \
        --remove-relation GEN-ass-TIS \
        --out data/processed/evaluation/graph_v4/benchmark_multiscale/ablation/graph_ablation_no_GEN-ass-TIS_v4.tsv
"""
import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-subgraph", required=True)
    ap.add_argument("--remove-relation", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df = pd.read_csv(args.base_subgraph, sep="\t", dtype=str)
    n_before = len(df)
    n_removed = int((df["relation"] == args.remove_relation).sum())
    if n_removed == 0:
        raise SystemExit(f"ERROR: la relacion '{args.remove_relation}' no aparece en {args.base_subgraph}")
    df_out = df[df["relation"] != args.remove_relation]
    df_out.to_csv(args.out, sep="\t", index=False)
    print(f"[ablation] relacion retirada: {args.remove_relation}")
    print(f"[ablation] aristas antes: {n_before:,} | retiradas: {n_removed:,} | despues: {len(df_out):,}")
    print(f"[ablation] -> {args.out}")


if __name__ == "__main__":
    main()
