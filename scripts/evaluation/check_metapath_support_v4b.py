#!/usr/bin/env python3
"""
Seccion 5 (TFM graph_v4): version 2 de la verificacion de soporte de
metapaths -- agrega, para cada metapath, percentiles de grado del hub
(p50/p75/p90/p95/p99) y una SIMULACION de cuanto quedaria el soporte
si se filtran los hubs por encima de cada percentil, antes de decidir
un umbral de poda concreto para GEN-TIS-GEN (que en la v1 mostro un
hub de grado 16.446 dominando el metapath -- ver metapath_support_v4.tsv).

Se aplica la misma simulacion a los 6 metapaths por consistencia/rigor,
no solo al que ya mostro el problema.

Uso:
    python3 check_metapath_support_v4b.py \
        --graph data/processed/integrated/graph_v4/graph_edges.tsv \
        --out data/processed/evaluation/graph_v4/metapath_support_v4_percentiles.tsv
"""
import argparse
from math import comb

import numpy as np
import pandas as pd

METAPATHS = [
    ("GEN-DIS-GEN", "GEN-ass-DIS", "node2"),
    ("GEN-TIS-GEN", "GEN-ass-TIS", "node2"),
    ("GEN-PWY-GEN", "GEN-ass-PWY", "node2"),
    ("GEN-CPD-GEN", "CPD-int-GEN", "node1"),
    ("GEN-CLL-GEN", "CLL-mut-GEN", "node1"),
    ("CPD-DIS-CPD", "CPD-trt-DIS", "node2"),
]

PERCENTILES = [50, 75, 90, 95, 99]


def log(msg):
    print(f"[metapath_support_v2] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    relations_needed = {rel for _, rel, _ in METAPATHS}
    log(f"Leyendo {args.graph} en chunks, filtrando relaciones: {sorted(relations_needed)}")

    chunks_by_rel = {rel: [] for rel in relations_needed}
    for chunk in pd.read_csv(args.graph, sep="\t", dtype=str,
                              usecols=["node1", "node2", "relation"],
                              chunksize=2_000_000):
        sub = chunk[chunk["relation"].isin(relations_needed)]
        for rel, df in sub.groupby("relation"):
            chunks_by_rel[rel].append(df)

    rows_out = []
    for name, rel, hub_col in METAPATHS:
        parts = chunks_by_rel.get(rel, [])
        if not parts:
            log(f"  {name}: relacion base '{rel}' no encontrada")
            continue
        df = pd.concat(parts, ignore_index=True)
        deg = df.groupby(hub_col).size()
        n_instancias_original = int(sum(comb(int(d), 2) for d in deg))

        pcts = {p: int(np.percentile(deg.values, p)) for p in PERCENTILES}
        log(f"  {name}: n_hubs={len(deg):,} percentiles_grado={pcts} "
            f"instancias_sin_filtrar={n_instancias_original:,}")

        row = {"metapath": name, "n_hubs": len(deg),
               "instancias_sin_filtrar": n_instancias_original}
        for p in PERCENTILES:
            umbral = pcts[p]
            deg_filtrado = deg[deg <= umbral]
            n_hubs_excluidos = len(deg) - len(deg_filtrado)
            n_instancias_filtradas = int(sum(comb(int(d), 2) for d in deg_filtrado))
            frac_retenida = (n_instancias_filtradas / n_instancias_original
                              if n_instancias_original > 0 else 0.0)
            row[f"p{p}_umbral_grado"] = umbral
            row[f"p{p}_hubs_excluidos"] = n_hubs_excluidos
            row[f"p{p}_instancias_resultantes"] = n_instancias_filtradas
            row[f"p{p}_frac_instancias_retenidas"] = round(frac_retenida, 4)
            log(f"    filtro p{p} (grado<= {umbral}): excluye {n_hubs_excluidos} hubs, "
                f"quedan {n_instancias_filtradas:,} instancias "
                f"({frac_retenida*100:.1f}% de las originales)")
        rows_out.append(row)

    pd.DataFrame(rows_out).to_csv(args.out, sep="\t", index=False)
    log(f"Resultado -> {args.out}")


if __name__ == "__main__":
    main()
