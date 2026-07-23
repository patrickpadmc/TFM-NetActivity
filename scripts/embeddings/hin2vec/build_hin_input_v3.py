#!/usr/bin/env python
"""
Filtra y prepara el edge set de entrada para HIN2Vec v3, a partir del
parquet cacheado con edge_class ya calculado (ver prepare_edge_class_v3.py).

4 variantes soportadas via --mode / --include-cex:
  auto + sin cex:      --mode auto
  auto + con cex:      --mode auto --include-cex
  metapath + sin cex:  --mode metapath --metapath-relations metapath_relations_v3.tsv
  metapath + con cex:  --mode metapath --metapath-relations metapath_relations_v3.tsv --include-cex

El submuestreo de GEN-cex-GEN usa la misma semilla y el mismo tamano
en auto+cex y metapath+cex para que ambas variantes sean comparables
entre si. Tamano por defecto: media de aristas del resto de relaciones
en graph_v3 (37 cargas de relacion, excluyendo cex) = 124,593.

Uso:
    python build_hin_input_v3.py --cached-parquet <path> --out <path> --mode auto [--include-cex]
    python build_hin_input_v3.py --cached-parquet <path> --out <path> --mode metapath \
        --metapath-relations metapath_relations_v3.tsv [--include-cex]
"""
import argparse
import time
import pandas as pd

CEX_CLASS = "GEN-cex-GEN"
DEFAULT_CEX_SAMPLE_SIZE = 124593  # media de aristas del resto de relaciones (excluyendo cex), graph_v3


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cached-parquet", required=True, help="Parquet de prepare_edge_class_v3.py")
    parser.add_argument("--out", required=True, help="Parquet de salida (edges listos para walker.py)")
    parser.add_argument("--mode", choices=["auto", "metapath"], required=True)
    parser.add_argument("--metapath-relations", default=None,
                         help="Whitelist de edge_class (requerido si --mode metapath)")
    parser.add_argument("--include-cex", action="store_true")
    parser.add_argument("--cex-sample-size", type=int, default=DEFAULT_CEX_SAMPLE_SIZE)
    parser.add_argument("--cex-sample-seed", type=int, default=42)
    args = parser.parse_args()

    if args.mode == "metapath" and not args.metapath_relations:
        parser.error("--metapath-relations es requerido si --mode metapath")

    log(f"Leyendo {args.cached_parquet}")
    df = pd.read_parquet(args.cached_parquet)
    log(f"{len(df)} aristas cacheadas")

    cex_mask = df["edge_class"] == CEX_CLASS
    non_cex = df[~cex_mask]
    cex = df[cex_mask]
    log(f"{len(cex)} aristas de cex, {len(non_cex)} del resto")

    if args.mode == "auto":
        base = non_cex
    else:
        with open(args.metapath_relations) as f:
            whitelist = set(line.strip() for line in f if line.strip())
        base = non_cex[non_cex["edge_class"].isin(whitelist)]
        log(f"Filtrado a metapath whitelist ({len(whitelist)} tipos): {len(base)} aristas")

    if args.include_cex:
        n = min(args.cex_sample_size, len(cex))
        cex_sample = cex.sample(n=n, random_state=args.cex_sample_seed)
        log(f"Incluyendo submuestra de cex: {len(cex_sample)} aristas (seed={args.cex_sample_seed})")
        final = pd.concat([base, cex_sample], ignore_index=True)
    else:
        final = base

    before = len(final)
    final = final.drop_duplicates(subset=["node1", "node2", "edge_class"])
    log(f"Deduplicado: {before} -> {len(final)} aristas")

    out_df = pd.DataFrame({
        "source_node": final["node1"],
        "source_class": final["node1_type"],
        "dest_node": final["node2"],
        "dest_class": final["node2_type"],
        "edge_class": final["edge_class"],
        "weight": 1.0,
    })
    out_df.to_parquet(args.out, index=False)
    log(f"Guardado {args.out}")
    log("Listo")


if __name__ == "__main__":
    main()
