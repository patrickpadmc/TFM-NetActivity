#!/usr/bin/env python3
"""
Seccion 10 (TFM graph_v4): extiende el catalogo unificado de embeddings
(Seccion 9) con las 4 condiciones nuevas entrenadas en esta seccion:
    gae / uniprot_fixed
    gae / uniprot_finetuned
    gat / uniprot_fixed
    gat / uniprot_finetuned

Mismo esquema que Seccion 9 (build_unified_embedding_catalog_v4.py):
    node_id | node_type | embedding_source | training_condition | seed | dimension | e_0001 ... e_n

Convierte SOLO los embeddings finales (75.128 nodos x 128-d), no los
vectores iniciales de UniProt proyectado (esos quedan como .npz aparte,
ya persistidos junto a los finales -- ver Seccion 10 seccion "Vectores
iniciales y finales"). Es un script separado del original de la
Seccion 9 (no se sobreescribe) para no perder la trazabilidad de que
epoca del proyecto genero cada archivo.

Uso:
    python3 build_unified_embedding_catalog_v4_seccion10.py \
        --subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --out-dir data/processed/analysis/embedding_catalog_v4
"""
import argparse
import os

import numpy as np
import pandas as pd


def log(msg):
    print(f"[catalog_s10] {msg}", flush=True)


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
    log(f"  {embedding_source}/{training_condition} (seed={seed}): {len(df):,} nodos x {dim} dim -> {out_path}"
        + (f"  [AVISO: {n_unknown} node_type no resueltos]" if n_unknown else ""))
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

    total = 0
    total += npz_to_unified(f"{EVAL}/gae_uniprot_fixed/gae_embeddings.npz", "gae", "uniprot_fixed", 42,
                             type_lookup, f"{args.out_dir}/gae_uniprot_fixed_seed42.tsv")
    total += npz_to_unified(f"{EVAL}/gae_uniprot_finetuned/gae_embeddings.npz", "gae", "uniprot_finetuned", 42,
                             type_lookup, f"{args.out_dir}/gae_uniprot_finetuned_seed42.tsv")
    total += npz_to_unified(f"{EVAL}/gat_uniprot_fixed/gat_embeddings.npz", "gat", "uniprot_fixed", 42,
                             type_lookup, f"{args.out_dir}/gat_uniprot_fixed_seed42.tsv")
    total += npz_to_unified(f"{EVAL}/gat_uniprot_finetuned/gat_embeddings.npz", "gat", "uniprot_finetuned", 42,
                             type_lookup, f"{args.out_dir}/gat_uniprot_finetuned_seed42.tsv")

    log(f"Total de filas escritas (4 condiciones nuevas): {total:,}")
    log("NOTA: los vectores iniciales de UniProt proyectado (*_initial_embeddings.npz, "
        "solo 21.111 nodos GEN cubiertos) NO se incluyen en este catalogo unificado -- "
        "quedan como .npz aparte junto a los finales, segun lo ya persistido en la Seccion 10.")
    log("Listo.")


if __name__ == "__main__":
    main()
