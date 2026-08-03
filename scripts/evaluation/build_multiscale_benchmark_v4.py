#!/usr/bin/env python3
"""
Seccion 5 (TFM graph_v4): construccion del subgrafo de benchmark
multi-escala (25/50/75/100%) para comparar los 5 metodos (Spectral,
Node2Vec, GAE, GAT-AE, HIN2Vec) de forma metodologicamente justa.

DECISION DE DISEÑO (confirmada por el usuario): el "100%" de la escala
NO es HIN2Vec's propio limite (que resulto ser barato incluso en el
grafo completo, ver job 308980: 6.24GB pico, ~2min) sino el techo real
impuesto por Spectral/GAE/GAT-AE, que en la Seccion 4 fallaron por
memoria en el grafo completo y solo funcionaron sobre el subgrafo comun
(~75.128 nodos, 3.339.427 aristas, ya construido y documentado en
build_common_subgraph_v4.py). Por lo tanto:
  100% = el subgrafo comun de la Seccion 4, TAL CUAL (reutilizado, no
         reconstruido desde cero -- ya tiene cobertura 100% verificada
         en los 6 tipos de nodo).
  75/50/25% = submuestras ESTRATIFICADAS ANIDADAS del 100% (no muestras
         independientes de graph_v4), preservando la proporcion relativa
         entre tipos de relacion que el 100% ya tiene, con la misma
         regla de dos pasos de la Seccion 4 para preservar conectividad
         (proteger aristas cuyo endpoint tiene grado 1, rescatar nodos
         que quedarian aislados tras el muestreo).

Uso:
    python3 build_multiscale_benchmark_v4.py \
        --base-subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --full-graph data/processed/integrated/graph_v4/graph_edges.tsv \
        --out-dir data/processed/evaluation/graph_v4/benchmark_multiscale \
        --seed 42
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

FRACTIONS = [0.25, 0.50, 0.75, 1.00]


def log(msg):
    print(f"[multiscale] {msg}", flush=True)


def stratified_sample_by_relation(df, fraction, seed):
    """Muestrea `fraction` de las aristas de CADA tipo de relacion por
    separado (estratificado), preservando la proporcion relativa entre
    relaciones que ya tiene df."""
    rng = np.random.default_rng(seed)
    parts = []
    for rel, group in df.groupby("relation"):
        n = max(1, int(round(len(group) * fraction)))
        idx = rng.choice(group.index.values, size=min(n, len(group)), replace=False)
        parts.append(group.loc[idx])
    return pd.concat(parts, ignore_index=True)


def protect_and_rescue(df_sampled, df_full, seed):
    """Regla de dos pasos (igual que Seccion 4):
    1) proteger aristas cuyo node1 o node2 tenga grado 1 en df_full
       (perderlas aislaria ese nodo de entrada).
    2) tras el muestreo, rescatar nodos que quedaron aislados agregando
       UNA arista candidata desde df_full que los reconecte."""
    deg_full = pd.concat([df_full["node1"], df_full["node2"]]).value_counts()
    deg1_nodes = set(deg_full[deg_full == 1].index)

    protect_mask = df_full["node1"].isin(deg1_nodes) | df_full["node2"].isin(deg1_nodes)
    protected = df_full[protect_mask]

    combined = pd.concat([df_sampled, protected], ignore_index=True).drop_duplicates(
        subset=["node1", "node2", "relation"])

    nodes_full = set(df_full["node1"]).union(df_full["node2"])
    nodes_sample = set(combined["node1"]).union(combined["node2"])
    stranded = nodes_full - nodes_sample
    if not stranded:
        return combined, 0

    candidate_mask = df_full["node1"].isin(stranded) | df_full["node2"].isin(stranded)
    candidates = df_full[candidate_mask]
    rescued_rows = []
    still_stranded = set(stranded)
    for _, row in candidates.iterrows():
        if row["node1"] in still_stranded or row["node2"] in still_stranded:
            rescued_rows.append(row)
            still_stranded.discard(row["node1"])
            still_stranded.discard(row["node2"])
        if not still_stranded:
            break
    if rescued_rows:
        combined = pd.concat([combined, pd.DataFrame(rescued_rows)], ignore_index=True) \
            .drop_duplicates(subset=["node1", "node2", "relation"])
    return combined, len(rescued_rows)


def coverage_report(df_sub, df_full):
    nodes_full = pd.concat([
        df_full[["node1", "node1_type"]].rename(columns={"node1": "node", "node1_type": "type"}),
        df_full[["node2", "node2_type"]].rename(columns={"node2": "node", "node2_type": "type"}),
    ]).drop_duplicates()
    nodes_sub = pd.concat([
        df_sub[["node1", "node1_type"]].rename(columns={"node1": "node", "node1_type": "type"}),
        df_sub[["node2", "node2_type"]].rename(columns={"node2": "node", "node2_type": "type"}),
    ]).drop_duplicates()
    cov_por_tipo = {}
    for t, grp in nodes_full.groupby("type"):
        n_full = len(grp)
        n_sub = nodes_sub[nodes_sub["type"] == t]["node"].nunique()
        cov_por_tipo[t] = {"n_full": int(n_full), "n_sub": int(n_sub),
                            "cobertura": round(n_sub / n_full, 4) if n_full else None}
    return cov_por_tipo, len(nodes_sub)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-subgraph", required=True,
                     help="subgrafo comun de la Seccion 4 (el '100%')")
    ap.add_argument("--full-graph", required=True,
                     help="graph_v4 completo, para calcular cobertura y para la regla de rescate")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    log(f"Cargando subgrafo base (100%) -> {args.base_subgraph}")
    df_base = pd.read_csv(args.base_subgraph, sep="\t", dtype=str)
    log(f"  {len(df_base):,} aristas, {df_base['relation'].nunique()} relaciones")

    log(f"Cargando graph_v4 completo (para cobertura/rescate) -> {args.full_graph}")
    full_chunks = []
    for chunk in pd.read_csv(args.full_graph, sep="\t", dtype=str, chunksize=5_000_000):
        full_chunks.append(chunk)
    df_full = pd.concat(full_chunks, ignore_index=True)
    log(f"  {len(df_full):,} aristas en graph_v4")

    manifest = {"seed": args.seed, "metodo": "muestreo estratificado por relacion, "
                "anidado dentro del subgrafo comun de la Seccion 4 (no independiente de graph_v4)",
                "fracciones": {}}

    for frac in FRACTIONS:
        tag = f"{int(frac * 100)}pct"
        log(f"--- Fraccion {tag} ---")
        if frac >= 1.0:
            df_final = df_base.copy()
            n_rescatados = 0
        else:
            df_sampled = stratified_sample_by_relation(df_base, frac, seed=args.seed)
            df_final, n_rescatados = protect_and_rescue(df_sampled, df_base, seed=args.seed)

        out_path = f"{args.out_dir}/graph_benchmark_{tag}_v4.tsv"
        df_final.to_csv(out_path, sep="\t", index=False)

        cov_por_tipo, n_nodos = coverage_report(df_final, df_full)
        log(f"  {len(df_final):,} aristas, {n_nodos:,} nodos, {n_rescatados} rescatados")
        for t, c in cov_por_tipo.items():
            log(f"    {t}: {c['n_sub']:,}/{c['n_full']:,} ({c['cobertura']*100:.1f}% de graph_v4)")

        manifest["fracciones"][tag] = {
            "archivo": out_path, "n_aristas": len(df_final), "n_nodos": n_nodos,
            "n_rescatados": n_rescatados, "cobertura_por_tipo": cov_por_tipo,
            "aristas_por_relacion": df_final["relation"].value_counts().to_dict(),
        }

    manifest_path = f"{args.out_dir}/benchmark_multiscale_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    log(f"Manifest -> {manifest_path}")


if __name__ == "__main__":
    main()
