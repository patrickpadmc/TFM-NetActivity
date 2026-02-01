from __future__ import annotations

from pathlib import Path
import json
import platform
import sys
import time

import pandas as pd
import networkx as nx


def main() -> None:
    t0 = time.time()

    root = Path(__file__).resolve().parents[2]  # repo root
    results_dir = root / "results"
    tables_dir = results_dir / "tables"
    figs_dir = results_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # --- 1) Create a tiny toy "edge list" dataset
    edges = pd.DataFrame(
        {
            "src": ["A", "A", "B", "C", "D", "E", "E"],
            "dst": ["B", "C", "D", "D", "E", "A", "C"],
            "weight": [1, 2, 1, 3, 1, 1, 2],
        }
    )

    raw_path = root / "data" / "raw" / "toy_edges.csv"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    edges.to_csv(raw_path, index=False)

    # --- 2) Load dataset back and build a weighted graph
    df = pd.read_csv(raw_path)
    G = nx.from_pandas_edgelist(df, source="src", target="dst", edge_attr="weight", create_using=nx.Graph)

    # --- 3) Compute simple graph metrics
    info = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": nx.density(G),
        "is_connected": nx.is_connected(G),
        "top_degree": sorted(G.degree, key=lambda x: x[1], reverse=True)[:5],
    }

    # --- 4) Save outputs
    summary_path = tables_dir / "smoke_test_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)

    # also save a small node table
    node_table = pd.DataFrame({"node": list(G.nodes()), "degree": [G.degree(n) for n in G.nodes()]})
    node_table_path = tables_dir / "smoke_test_nodes.csv"
    node_table.to_csv(node_table_path, index=False)

    meta = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "elapsed_s": round(time.time() - t0, 3),
        "repo_root": str(root),
    }
    meta_path = tables_dir / "smoke_test_meta.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print("SMOKE TEST OK")
    print(f"- wrote: {raw_path}")
    print(f"- wrote: {summary_path}")
    print(f"- wrote: {node_table_path}")
    print(f"- wrote: {meta_path}")


if __name__ == "__main__":
    main()