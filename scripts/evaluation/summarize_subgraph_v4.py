#!/usr/bin/env python3
"""
Seccion 7 (TFM graph_v4): agrega el subgrafo comun por combinacion
base de datos + relacion (mismo esquema que la Tabla 1 de Seccion 1),
para poder construir subgraph_summary.csv con las mismas columnas.
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subgraph", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.subgraph, sep="\t", dtype=str)
    print(f"[summarize] {len(df):,} aristas totales en el subgrafo")

    rows = []
    for (db, rel), g in df.groupby(["database", "relation"]):
        n_aristas = len(g)
        nodos_unicos = pd.concat([g["node1"], g["node2"]]).nunique()
        rows.append({"database": db, "relation": rel, "n_aristas": n_aristas, "nodos_unicos": nodos_unicos})

    out_df = pd.DataFrame(rows).sort_values(["database", "relation"]).reset_index(drop=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[summarize] {len(out_df)} combinaciones fuente+relacion -> {args.out}")
    print(f"[summarize] suma aristas: {out_df['n_aristas'].sum():,}")

    nodos_unicos_totales = pd.concat([df["node1"], df["node2"]]).nunique()
    print(f"[summarize] nodos unicos totales del subgrafo: {nodos_unicos_totales:,}")


if __name__ == "__main__":
    main()
