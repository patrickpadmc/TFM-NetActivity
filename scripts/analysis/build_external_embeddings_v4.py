#!/usr/bin/env python3
"""
Seccion 8 (TFM graph_v4): integracion y auditoria de embeddings externos
(UniProt ProtT5, GenePT-Model3) para los nodos GEN del subgrafo comun.

Estrategia de seleccion canonica (documentada, nunca silenciosa/aleatoria):
  1. Candidatos UniProt de un ENSG = filas de ens2uniprot.tsv (Bioteque).
  2. Se filtran por 'reviewed' (presentes en human_reviewed.tsv, Swiss-Prot
     humano) -- criterio biologico: entrada revisada = anotacion curada
     manualmente, preferible a TrEMBL no revisada.
  3. Si hay exactamente 1 candidato revisado -> canonico.
  4. Si hay >1 candidato revisado (ambiguedad real) -> desempate
     deterministico: (a) el que tiene embedding disponible en el .h5;
     si sigue empatado, (b) orden alfabetico de accesion. Nunca al azar.
  5. Si hay 0 revisados pero >=1 no revisado -> fallback documentado,
     estado_mapeo = 'mapeado_no_revisado' (no se descarta, se marca).
     Si hay mas de un no revisado sin forma de desempatar (ninguno con
     embedding), se excluye por ambiguedad no resoluble.
  6. Si no hay ningun candidato -> sin_mapeo (motivo: 'sin entrada en
     ens2uniprot').

Bridge ENSG -> simbolo genico: no existe un mapeo directo ENSG->simbolo
en los recursos de Bioteque disponibles, asi que se deriva vía la
accesion UniProt canonica ya seleccionada: ENSG -> UniProt -> simbolo,
usando gname2uniprot.tsv invertido (simbolo->UniProt, se invierte a
UniProt->simbolo). Esto se documenta explicitamente porque es una
decision metodologica, no un mapeo directo de una fuente unica.

GenePT-Model3 se busca por simbolo (mayusculas) directamente en el
diccionario descargado de Zenodo.

Uso:
    python3 build_external_embeddings_v4.py \
        --subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --ens2uniprot external/../GEN/ens2uniprot.tsv \
        --reviewed external/../GEN/human_reviewed.tsv \
        --gname2uniprot external/../GEN/gname2uniprot.tsv \
        --uniprot-h5 data/raw/databases/uniprot_embeddings/UP000005640_9606/per-protein.h5 \
        --genept-pickle data/raw/databases/genept/GenePT_gene_protein_embedding_model_3_text.pickle \
        --out-dir data/processed/analysis/external_embeddings_v4
"""
import argparse
import json
import os

import h5py
import numpy as np
import pandas as pd
import pickle


def log(msg):
    print(f"[ext_embeddings] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subgraph", required=True)
    ap.add_argument("--ens2uniprot", required=True)
    ap.add_argument("--reviewed", required=True)
    ap.add_argument("--gname2uniprot", required=True)
    ap.add_argument("--uniprot-h5", required=True)
    ap.add_argument("--genept-pickle", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Genes GEN unicos del subgrafo
    log(f"Cargando subgrafo -> {args.subgraph}")
    edges = pd.read_csv(args.subgraph, sep="\t", dtype=str)
    gen_nodes = set(edges.loc[edges["node1_type"] == "GEN", "node1"]) | \
                set(edges.loc[edges["node2_type"] == "GEN", "node2"])
    gen_nodes = sorted(gen_nodes)
    log(f"GEN unicos en el subgrafo: {len(gen_nodes):,}")

    # 2. Recursos de mapeo
    log("Cargando ens2uniprot, human_reviewed, gname2uniprot...")
    ens2up = pd.read_csv(args.ens2uniprot, sep=r"\s+", header=None,
                          names=["ensg", "uniprot"], engine="python")
    ens2up_map = ens2up.groupby("ensg")["uniprot"].apply(list).to_dict()

    reviewed_df = pd.read_csv(args.reviewed, sep="\t")
    reviewed_set = set(reviewed_df["Entry"])
    log(f"  accesiones revisadas (Swiss-Prot humano): {len(reviewed_set):,}")

    gname2up = pd.read_csv(args.gname2uniprot, sep="\t", header=None,
                            names=["symbol", "uniprot"])
    up2symbol = {}
    n_up_multi_symbol = 0
    for up, grp in gname2up.groupby("uniprot"):
        symbols = sorted(grp["symbol"].unique())
        if len(symbols) > 1:
            n_up_multi_symbol += 1
        up2symbol[up] = symbols[0]
    log(f"  accesiones con >1 simbolo en gname2uniprot (se tomo el primero alfabetico): {n_up_multi_symbol}")

    with h5py.File(args.uniprot_h5, "r") as f:
        uniprot_emb_keys = set(f.keys())
    log(f"  accesiones con embedding UniProt disponible: {len(uniprot_emb_keys):,}")

    log(f"Cargando GenePT-Model3 pickle (puede tardar, ~900MB)...")
    with open(args.genept_pickle, "rb") as f:
        genept_dict = pickle.load(f)
    genept_keys_upper = {k.upper() for k in genept_dict.keys()}
    log(f"  simbolos con embedding GenePT-Model3: {len(genept_dict):,}")

    # 3. Seleccion canonica por ENSG
    log("Aplicando estrategia de seleccion canonica...")
    rows = []
    for ensg in gen_nodes:
        candidatos = ens2up_map.get(ensg, [])
        row = {"graph_node_id": ensg, "identificador_gen_original": ensg,
               "simbolo_genico": None, "uniprot_accession": None,
               "estado_mapeo": None, "motivo_exclusion": None}

        if not candidatos:
            row["estado_mapeo"] = "sin_mapeo"
            row["motivo_exclusion"] = "sin entrada en ens2uniprot (Bioteque)"
            rows.append(row)
            continue

        revisados = [c for c in candidatos if c in reviewed_set]
        no_revisados = [c for c in candidatos if c not in reviewed_set]

        elegido = None
        if len(revisados) == 1:
            elegido = revisados[0]
            row["estado_mapeo"] = "mapeado_revisado"
        elif len(revisados) > 1:
            con_embedding = [c for c in revisados if c in uniprot_emb_keys]
            if len(con_embedding) >= 1:
                elegido = sorted(con_embedding)[0]
            else:
                elegido = sorted(revisados)[0]
            row["estado_mapeo"] = "mapeado_revisado_desempate"
            row["motivo_exclusion"] = f"multiples revisados ({len(revisados)}), desempate por embedding disponible + orden alfabetico"
        elif len(no_revisados) == 1:
            elegido = no_revisados[0]
            row["estado_mapeo"] = "mapeado_no_revisado"
            row["motivo_exclusion"] = "sin entrada Swiss-Prot revisada; se usa la unica accesion TrEMBL disponible"
        elif len(no_revisados) > 1:
            con_embedding = [c for c in no_revisados if c in uniprot_emb_keys]
            if len(con_embedding) == 1:
                elegido = con_embedding[0]
                row["estado_mapeo"] = "mapeado_no_revisado_desempate"
                row["motivo_exclusion"] = f"multiples no revisados ({len(no_revisados)}), desempate por unico con embedding disponible"
            else:
                row["estado_mapeo"] = "excluido"
                row["motivo_exclusion"] = f"multiples accesiones no revisadas ({len(no_revisados)}) sin criterio de desempate resoluble (0 o >1 con embedding)"

        if elegido is not None:
            row["uniprot_accession"] = elegido
            row["simbolo_genico"] = up2symbol.get(elegido)

        rows.append(row)

    map_df = pd.DataFrame(rows)

    # 4. Disponibilidad de embedding UniProt real (accesion elegida está en el .h5)
    map_df["uniprot_has_embedding"] = map_df["uniprot_accession"].isin(uniprot_emb_keys)

    # 5. Disponibilidad GenePT-Model3 (via simbolo, case-insensitive)
    map_df["genept_has_embedding"] = map_df["simbolo_genico"].apply(
        lambda s: (s is not None) and (str(s).upper() in genept_keys_upper)
    )

    map_path = f"{args.out_dir}/gen_id_mapping_table.tsv"
    map_df.to_csv(map_path, sep="\t", index=False)
    log(f"Tabla de mapeo -> {map_path}")

    # 6. Resumen de cobertura global
    n_total = len(map_df)
    n_uniprot_ok = int(map_df["uniprot_has_embedding"].sum())
    n_genept_ok = int(map_df["genept_has_embedding"].sum())

    summary = {
        "n_genes_subgrafo": n_total,
        "uniprot": {
            "n_mapeados_con_embedding": n_uniprot_ok,
            "pct_mapeados_con_embedding": round(100 * n_uniprot_ok / n_total, 2),
            "n_sin_embedding": n_total - n_uniprot_ok,
            "pct_sin_embedding": round(100 * (n_total - n_uniprot_ok) / n_total, 2),
        },
        "genept_model3": {
            "n_mapeados_con_embedding": n_genept_ok,
            "pct_mapeados_con_embedding": round(100 * n_genept_ok / n_total, 2),
            "n_sin_embedding": n_total - n_genept_ok,
            "pct_sin_embedding": round(100 * (n_total - n_genept_ok) / n_total, 2),
        },
        "estado_mapeo_breakdown": map_df["estado_mapeo"].value_counts().to_dict(),
    }
    with open(f"{args.out_dir}/coverage_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    log(f"Cobertura global UniProt: {n_uniprot_ok}/{n_total} ({summary['uniprot']['pct_mapeados_con_embedding']}%)")
    log(f"Cobertura global GenePT-Model3: {n_genept_ok}/{n_total} ({summary['genept_model3']['pct_mapeados_con_embedding']}%)")

    # 7. Cobertura por relacion que incluya GEN
    log("Calculando cobertura por relacion...")
    ensg2uniprot_ok = dict(zip(map_df["graph_node_id"], map_df["uniprot_has_embedding"]))
    ensg2genept_ok = dict(zip(map_df["graph_node_id"], map_df["genept_has_embedding"]))

    rel_rows = []
    gen_edges = edges[(edges["node1_type"] == "GEN") | (edges["node2_type"] == "GEN")]
    for rel, g in gen_edges.groupby("relation"):
        gen_ids_in_rel = set()
        gen_ids_in_rel |= set(g.loc[g["node1_type"] == "GEN", "node1"])
        gen_ids_in_rel |= set(g.loc[g["node2_type"] == "GEN", "node2"])
        n_gen = len(gen_ids_in_rel)
        n_up = sum(1 for gid in gen_ids_in_rel if ensg2uniprot_ok.get(gid, False))
        n_gp = sum(1 for gid in gen_ids_in_rel if ensg2genept_ok.get(gid, False))
        rel_rows.append({
            "relation": rel, "n_aristas": len(g), "n_genes_unicos": n_gen,
            "n_uniprot_cobertura": n_up, "pct_uniprot_cobertura": round(100 * n_up / n_gen, 2) if n_gen else None,
            "n_genept_cobertura": n_gp, "pct_genept_cobertura": round(100 * n_gp / n_gen, 2) if n_gen else None,
        })
    rel_df = pd.DataFrame(rel_rows).sort_values("n_aristas", ascending=False)
    rel_path = f"{args.out_dir}/coverage_by_relation.tsv"
    rel_df.to_csv(rel_path, sep="\t", index=False)
    log(f"Cobertura por relacion -> {rel_path}")

    # 8. Exportar embeddings normalizados (solo genes con embedding real)
    log("Exportando embeddings UniProt normalizados...")
    up_rows = map_df[map_df["uniprot_has_embedding"]]
    with h5py.File(args.uniprot_h5, "r") as f:
        vectors = []
        node_ids = []
        for _, r in up_rows.iterrows():
            vectors.append(np.asarray(f[r["uniprot_accession"]][:], dtype=np.float32))
            node_ids.append(r["graph_node_id"])
    dim_up = len(vectors[0]) if vectors else 0
    up_out = pd.DataFrame(vectors, columns=[f"e_{i+1:04d}" for i in range(dim_up)])
    up_out.insert(0, "node_id", node_ids)
    up_out.insert(1, "node_type", "GEN")
    up_out.insert(2, "embedding_source", "UniProt_ProtT5")
    up_out.insert(3, "dimension", dim_up)
    up_out.to_csv(f"{args.out_dir}/uniprot_embeddings_subgraph.tsv", sep="\t", index=False)
    log(f"  {len(up_out)} genes x {dim_up} dim -> uniprot_embeddings_subgraph.tsv")

    log("Exportando embeddings GenePT-Model3 normalizados...")
    gp_rows = map_df[map_df["genept_has_embedding"]]
    genept_upper_lookup = {k.upper(): k for k in genept_dict.keys()}
    vectors = []
    node_ids = []
    for _, r in gp_rows.iterrows():
        real_key = genept_upper_lookup[str(r["simbolo_genico"]).upper()]
        vectors.append(np.asarray(genept_dict[real_key], dtype=np.float32))
        node_ids.append(r["graph_node_id"])
    dim_gp = len(vectors[0]) if vectors else 0
    gp_out = pd.DataFrame(vectors, columns=[f"e_{i+1:04d}" for i in range(dim_gp)])
    gp_out.insert(0, "node_id", node_ids)
    gp_out.insert(1, "node_type", "GEN")
    gp_out.insert(2, "embedding_source", "GenePT_Model3")
    gp_out.insert(3, "dimension", dim_gp)
    gp_out.to_csv(f"{args.out_dir}/genept_embeddings_subgraph.tsv", sep="\t", index=False)
    log(f"  {len(gp_out)} genes x {dim_gp} dim -> genept_embeddings_subgraph.tsv")

    # 9. Informe de exclusion
    excl_df = map_df[~map_df["uniprot_has_embedding"] | ~map_df["genept_has_embedding"]][
        ["graph_node_id", "simbolo_genico", "uniprot_accession", "estado_mapeo",
         "motivo_exclusion", "uniprot_has_embedding", "genept_has_embedding"]
    ]
    excl_path = f"{args.out_dir}/genes_excluidos_report.tsv"
    excl_df.to_csv(excl_path, sep="\t", index=False)
    log(f"Genes con al menos una fuente sin embedding -> {excl_path} ({len(excl_df)} filas)")

    log("Listo.")


if __name__ == "__main__":
    main()
