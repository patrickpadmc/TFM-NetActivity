#!/usr/bin/env python3
"""
Seccion 9 (TFM graph_v4): convierte los .npz/.tsv ya existentes de las
5 condiciones estructurales + 2 crudas externas al esquema de tabla
unificado exigido por el protocolo:

    node_id | node_type | embedding_source | training_condition | seed | dimension | e_0001 ... e_n

Incluye tambien las replicas de semilla (43/44) que ya existan para
GAE y HIN2Vec (Seccion 6). Node2Vec seed 43/44 NO se incluye porque
esas corridas son anteriores al parche de persistencia (ver informe).

Uso:
    python3 build_unified_embedding_catalog_v4.py \
        --subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --out-dir data/processed/analysis/embedding_catalog_v4
"""
import argparse
import os

import numpy as np
import pandas as pd


def log(msg):
    print(f"[catalog] {msg}", flush=True)


def node_type_lookup(subgraph_path):
    df = pd.read_csv(subgraph_path, sep="\t", dtype=str,
                      usecols=["node1", "node2", "node1_type", "node2_type"])
    lookup = {}
    for n, t in zip(df["node1"], df["node1_type"]):
        lookup[n] = t
    for n, t in zip(df["node2"], df["node2_type"]):
        lookup[n] = t
    return lookup


def npz_to_unified(npz_path, embedding_source, training_condition, seed, type_lookup, out_path):
    data = np.load(npz_path, allow_pickle=True)
    node_ids = data["node_ids"]
    embeddings = data["embeddings"]
    dim = embeddings.shape[1]
    df = pd.DataFrame(embeddings, columns=[f"e_{i+1:04d}" for i in range(dim)])
    df.insert(0, "node_id", node_ids)
    df.insert(1, "node_type", [type_lookup.get(n, "DESCONOCIDO") for n in node_ids])
    df.insert(2, "embedding_source", embedding_source)
    df.insert(3, "training_condition", training_condition)
    df.insert(4, "seed", seed)
    df.insert(5, "dimension", dim)
    df.to_csv(out_path, sep="\t", index=False)
    n_unknown = (df["node_type"] == "DESCONOCIDO").sum()
    log(f"  {embedding_source} (seed={seed}): {len(df):,} nodos x {dim} dim -> {out_path}"
        + (f"  [AVISO: {n_unknown} node_type no resueltos]" if n_unknown else ""))
    return len(df)


def raw_tsv_to_unified(tsv_path, training_condition, seed, out_path):
    df = pd.read_csv(tsv_path, sep="\t")
    df.insert(3, "training_condition", training_condition)
    df.insert(4, "seed", seed)
    df.to_csv(out_path, sep="\t", index=False)
    log(f"  {df['embedding_source'].iloc[0]} (crudo): {len(df):,} nodos -> {out_path}")
    return len(df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subgraph", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    log("Construyendo lookup node_id -> node_type desde el subgrafo...")
    type_lookup = node_type_lookup(args.subgraph)
    log(f"  {len(type_lookup):,} nodos unicos indexados")

    EVAL = "data/processed/evaluation/graph_v4"
    BENCH = f"{EVAL}/benchmark_multiscale"
    EXT = "data/processed/analysis/external_embeddings_v4"

    total = 0
    total += npz_to_unified(f"{EVAL}/embeddings/spectral_v4_subgraph.npz", "spectral", "structural", 42,
                             type_lookup, f"{args.out_dir}/spectral_structural_seed42.tsv")
    total += npz_to_unified(f"{BENCH}/node2vec/100pct/node2vec_embeddings.npz", "node2vec", "structural", 42,
                             type_lookup, f"{args.out_dir}/node2vec_structural_seed42.tsv")
    total += npz_to_unified(f"{BENCH}/hin2vec/100pct/hin2vec_embeddings.npz", "hin2vec", "structural", 42,
                             type_lookup, f"{args.out_dir}/hin2vec_structural_seed42.tsv")
    total += npz_to_unified(f"{EVAL}/gae/gae_embeddings.npz", "gae", "structural", 42,
                             type_lookup, f"{args.out_dir}/gae_structural_seed42.tsv")
    total += npz_to_unified(f"{EVAL}/gat/gat_embeddings.npz", "gat", "structural", 42,
                             type_lookup, f"{args.out_dir}/gat_structural_seed42.tsv")

    for seed in (43, 44):
        gae_path = f"{BENCH}/gae/100pct_seed{seed}/gae_embeddings.npz"
        if os.path.exists(gae_path):
            total += npz_to_unified(gae_path, "gae", "structural", seed,
                                     type_lookup, f"{args.out_dir}/gae_structural_seed{seed}.tsv")
        else:
            log(f"  AVISO: no existe {gae_path}, se omite gae seed={seed}")

        hin_path = f"{BENCH}/hin2vec/100pct_seed{seed}/hin2vec_embeddings.npz"
        if os.path.exists(hin_path):
            total += npz_to_unified(hin_path, "hin2vec", "structural", seed,
                                     type_lookup, f"{args.out_dir}/hin2vec_structural_seed{seed}.tsv")
        else:
            log(f"  AVISO: no existe {hin_path}, se omite hin2vec seed={seed}")

    log("AVISO: node2vec seed=43/44 NO se incluye -- esas corridas de la Seccion 6 son anteriores "
        "al parche de persistencia de embeddings; solo se guardaron metricas, no vectores.")

    total += raw_tsv_to_unified(f"{EXT}/uniprot_embeddings_subgraph.tsv", "raw_external", "N/A",
                                 f"{args.out_dir}/uniprot_raw_external.tsv")
    total += raw_tsv_to_unified(f"{EXT}/genept_embeddings_subgraph.tsv", "raw_external", "N/A",
                                 f"{args.out_dir}/genept_raw_external.tsv")

    log(f"Total de filas escritas en el catalogo unificado: {total:,}")
    log("Listo.")


if __name__ == "__main__":
    main()
