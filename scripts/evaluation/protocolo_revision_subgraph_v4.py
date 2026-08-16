#!/usr/bin/env python3
import pandas as pd

p = "data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv"
df = pd.read_csv(p, sep="\t", dtype=str)

nodes = set(df.node1) | set(df.node2)
degree = pd.concat([df.node1, df.node2]).value_counts()

print(f"\nARISTAS={len(df)}")
print(f"NODOS={len(nodes)}")
print(f"RELACIONES={df.relation.nunique()}")
print(f"NODOS_GRADO_0={sum(n not in degree.index for n in nodes)}")
print(f"NODOS_GRADO_1={(degree == 1).sum()}")

print("\nAristas por relación:")
print(df.groupby("relation").size().sort_values(ascending=False).to_string())
