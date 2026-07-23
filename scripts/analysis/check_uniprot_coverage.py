#!/usr/bin/env python
"""
Calcula la cobertura de mapeo ENSG -> UniProt y la disponibilidad de
embeddings UniProt (per-protein.h5) para los genes presentes en graph_v2.

Uso:
    python check_uniprot_coverage.py \
        --graph-edges data/processed/integrated/graph_v2/graph_edges.tsv \
        --ens2uniprot external/bioteque/metadata/mappings/GEN/ens2uniprot.tsv \
        --embeddings-h5 data/raw/databases/uniprot_embeddings/UP000005640_9606/per-protein.h5 \
        --out-dir data/processed/analysis/uniprot_coverage

Importante: ejecutar en nodo de computo (srun/sbatch), no en login node,
porque usa pandas.
"""

import argparse
from pathlib import Path

import h5py
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--graph-edges", required=True,
                    help="TSV con columnas node1,node2,node1_type,node2_type,relation,database")
    p.add_argument("--ens2uniprot", required=True,
                    help="TSV de Bioteque con columnas ensg, uniprot (sin header)")
    p.add_argument("--embeddings-h5", required=True,
                    help="Archivo per-protein.h5 de UniProt (keys = accessions)")
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def extract_gen_nodes(graph_edges_path):
    edges = pd.read_csv(graph_edges_path, sep="\t")
    gen_nodes = set()
    gen_nodes |= set(edges.loc[edges["node1_type"] == "GEN", "node1"])
    gen_nodes |= set(edges.loc[edges["node2_type"] == "GEN", "node2"])
    return sorted(gen_nodes)


def load_mapping(mapping_path):
    # sep="\s+" tolera tanto tab como espacio simple, por si el delimitador
    # real difiere de lo esperado
    mapping = pd.read_csv(
        mapping_path, sep=r"\s+", header=None, names=["ensg", "uniprot"], engine="python"
    )
    return mapping


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Genes (ENSG) presentes en el grafo v2
    gen_nodes = extract_gen_nodes(args.graph_edges)
    n_ensg_total = len(gen_nodes)
    print(f"GEN unicos en graph_v2: {n_ensg_total}")

    # 2. Mapeo ENSG -> UniProt (puede ser 1-a-muchos)
    mapping = load_mapping(args.ens2uniprot)
    map_dict = mapping.groupby("ensg")["uniprot"].apply(list).to_dict()

    rows = []
    for ensg in gen_nodes:
        uniprots = map_dict.get(ensg, [])
        if uniprots:
            for up in uniprots:
                rows.append({"ensg": ensg, "uniprot": up})
        else:
            rows.append({"ensg": ensg, "uniprot": None})

    mapped_df = pd.DataFrame(rows)

    mapped_mask = mapped_df["uniprot"].notna()
    n_rows_post_mapeo = int(mapped_mask.sum())
    n_ensg_unicos_post_mapeo = mapped_df.loc[mapped_mask, "ensg"].nunique()
    n_uniprot_unicos_post_mapeo = mapped_df.loc[mapped_mask, "uniprot"].nunique()

    # 3. Disponibilidad en el archivo de embeddings
    with h5py.File(args.embeddings_h5, "r") as f:
        embedding_keys = set(f.keys())

    mapped_df["has_embedding"] = mapped_df["uniprot"].isin(embedding_keys)
    emb_df = mapped_df[mapped_df["has_embedding"]]

    n_rows_con_embedding = len(emb_df)
    n_ensg_unicos_con_embedding = emb_df["ensg"].nunique()
    n_uniprot_unicos_con_embedding = emb_df["uniprot"].nunique()

    # 4. Resumen
    summary_rows = [
        ("GEN en grafo v2 (ENSG)", n_ensg_total, n_ensg_total, "-"),
        ("Con mapeo a UniProt", n_rows_post_mapeo, n_ensg_unicos_post_mapeo, n_uniprot_unicos_post_mapeo),
        ("Con embedding disponible en .h5", n_rows_con_embedding, n_ensg_unicos_con_embedding, n_uniprot_unicos_con_embedding),
    ]

    summary_df = pd.DataFrame(
        summary_rows, columns=["etapa", "filas_totales", "ensg_unicos", "uniprot_unicos"]
    )

    # Guardar outputs
    summary_path = out_dir / "uniprot_embedding_coverage_summary.tsv"
    full_mapping_path = out_dir / "ensg_uniprot_mapping_full.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    mapped_df.to_csv(full_mapping_path, sep="\t", index=False)

    # Tabla markdown en stdout
    print()
    print("|Etapa|Filas totales|ENSG unicos|UniProt accessions unicos|")
    print("|---|---|---|---|")
    for etapa, filas, ensg_u, up_u in summary_rows:
        print(f"|{etapa}|{filas}|{ensg_u}|{up_u}|")

    print()
    print(f"Resumen guardado en: {summary_path}")
    print(f"Mapeo completo guardado en: {full_mapping_path}")

    # Aviso si hay multi-mapeo ENSG -> multiples UniProt
    n_multi = mapped_df.groupby("ensg")["uniprot"].apply(lambda x: x.notna().sum() > 1).sum()
    if n_multi > 0:
        print(f"\nAVISO: {n_multi} genes ENSG mapean a mas de un UniProt accession (isoformas).")


if __name__ == "__main__":
    main()
