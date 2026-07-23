#!/usr/bin/env python3
"""
count_nodes_edges_by_type.py
Calcula, para el grafo (graph_edges.tsv):
  - Nº de Nodos: nodos unicos por tipo (DIS, GEN, TIS, PWY, CPD, CLL)
  - Matriz Nº de Aristas <tipo1> x <tipo2>: aristas entre cada par de
    tipos (incluida la diagonal, p.ej. GEN-GEN). Una arista T1-T2 cuenta
    una vez en la celda [T1][T2] y, por simetria, es la misma celda que
    [T2][T1] (la matriz es simetrica porque node1/node2 son solo el
    orden en que se guardo la arista, no una direccion logica).

Input:  graph_edges.tsv (columnas node1, node2, node1_type, node2_type,
        relation, database)
Output: TSV con la tabla + impresion en pantalla en formato markdown

Uso:
  python3 count_nodes_edges_by_type.py
  python3 count_nodes_edges_by_type.py --graph-edges <ruta> --out <ruta>
"""
import argparse
import os
import pandas as pd

NODE_TYPES = ["DIS", "GEN", "TIS", "PWY", "CPD", "CLL"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--graph-edges",
        default=os.path.expanduser(
            "~/work/TFM-NetActivity/data/processed/integrated/graph_v2/graph_edges.tsv"
        ),
    )
    p.add_argument(
        "--out",
        default=os.path.expanduser(
            "~/work/TFM-NetActivity/data/processed/integrated/graph_v2/node_edge_stats.tsv"
        ),
    )
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.graph_edges, sep="\t", dtype=str)

    # Nº de nodos unicos por tipo
    n_nodes = {}
    for t in NODE_TYPES:
        nodes_as_n1 = set(df.loc[df["node1_type"] == t, "node1"])
        nodes_as_n2 = set(df.loc[df["node2_type"] == t, "node2"])
        n_nodes[t] = len(nodes_as_n1 | nodes_as_n2)

    # Matriz tipo x tipo de aristas (simetrica, incluida la diagonal)
    matrix = {t1: {} for t1 in NODE_TYPES}
    for t1 in NODE_TYPES:
        for t2 in NODE_TYPES:
            if t1 == t2:
                mask = (df["node1_type"] == t1) & (df["node2_type"] == t1)
            else:
                mask = (
                    ((df["node1_type"] == t1) & (df["node2_type"] == t2))
                    | ((df["node1_type"] == t2) & (df["node2_type"] == t1))
                )
            matrix[t1][t2] = int(mask.sum())

    # Construir filas de salida
    rows = []
    for t1 in NODE_TYPES:
        row = {"Tipo Nodo": t1, "Nº Nodos": n_nodes[t1]}
        for t2 in NODE_TYPES:
            row[f"Nº Aristas {t2}"] = matrix[t1][t2]
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, sep="\t", index=False)

    # Tabla markdown para pegar directamente en la memoria/informe
    header = "|Tipo Nodo|Nº Nodos|" + "|".join(f"Nº Aristas {t}" for t in NODE_TYPES) + "|"
    sep = "|---------|--------|" + "|".join("-" * len(f"Nº Aristas {t}") for t in NODE_TYPES) + "|"
    print()
    print(header)
    print(sep)
    for row in rows:
        cells = [row["Tipo Nodo"], str(row["Nº Nodos"])] + [
            str(row[f"Nº Aristas {t}"]) for t in NODE_TYPES
        ]
        print("|" + "|".join(cells) + "|")

    print(f"\nTSV guardado en: {args.out}")


if __name__ == "__main__":
    main()
