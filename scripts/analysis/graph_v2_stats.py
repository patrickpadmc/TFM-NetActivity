#!/usr/bin/env python
"""
Calcula estadisticas basicas (nodos y aristas) de graph_v2.

Uso:
    python graph_v2_stats.py \
        --graph-edges data/processed/integrated/graph_v2/graph_edges.tsv \
        --out-dir data/processed/analysis/graph_v2_stats

Ejecutar en nodo de computo (srun/sbatch), no en login node, por uso de pandas.
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--graph-edges", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(args.graph_edges, sep="\t")

    n_edges_filas = len(edges)
    n_edges_unicas = edges.drop_duplicates(subset=["node1", "node2", "relation"]).shape[0]

    # Nodos: combinar (node1,node1_type) y (node2,node2_type)
    nodes1 = edges[["node1", "node1_type"]].rename(
        columns={"node1": "node_id", "node1_type": "node_type"}
    )
    nodes2 = edges[["node2", "node2_type"]].rename(
        columns={"node2": "node_id", "node2_type": "node_type"}
    )
    all_nodes = pd.concat([nodes1, nodes2], ignore_index=True).drop_duplicates()

    n_nodes_total = all_nodes["node_id"].nunique()
    nodes_by_type = all_nodes.groupby("node_type")["node_id"].nunique().sort_values(ascending=False)

    summary_path = out_dir / "graph_v2_stats_summary.tsv"
    nodes_by_type_path = out_dir / "graph_v2_nodes_by_type.tsv"

    with open(summary_path, "w") as f:
        f.write(f"filas_totales_aristas\t{n_edges_filas}\n")
        f.write(f"aristas_unicas\t{n_edges_unicas}\n")
        f.write(f"nodos_unicos_totales\t{n_nodes_total}\n")

    nodes_by_type.to_csv(nodes_by_type_path, sep="\t", header=["n_nodos"])

    print(f"Filas totales (aristas): {n_edges_filas}")
    print(f"Aristas unicas (node1,node2,relation): {n_edges_unicas}")
    print(f"Nodos unicos totales: {n_nodes_total}")
    print()
    print("Nodos por tipo:")
    print(nodes_by_type.to_string())

    if n_edges_filas != n_edges_unicas:
        n_dup = n_edges_filas - n_edges_unicas
        print(
            f"\nAVISO: {n_dup} filas duplicadas en (node1,node2,relation) "
            f"- revisar si vienen de multiples databases reportando la misma relacion."
        )


if __name__ == "__main__":
    main()
