#!/usr/bin/env python3
"""
Seccion 4/5 (TFM graph_v4): subgrafo comun.

Motivado por la Seccion 4: Spectral Embedding sobre el grafo de
entrenamiento completo (110,324,395 aristas) fallo por falta de memoria
en la factorizacion LU dispersa necesaria para el modo shift-invert de
eigsh (job 308475, MemoryError con 450GB asignados). Se construye aqui
el subgrafo comun (definido en la Seccion 5 del protocolo) para
reintentar Spectral, y para los demas metodos de la Seccion 4 si
tambien resultan inviables sobre el grafo completo.

Requisitos del subgrafo (Seccion 5):
  - Conservar los 6 tipos de nodo cuando sea posible.
  - Conservar las relaciones primarias seleccionadas (Seccion 2).
  - Mantener de forma aproximada la proporcion de nodos/aristas por tipo.
  - Ser viable para HIN2Vec.
  - Muestreo estratificado, nunca seleccion arbitraria de nodos.

Diseno de muestreo (confirmado con el usuario, ver docs/graph_v4/):
  - Relaciones primarias (7, ya train-only): 100% -- son la senal de
    evaluacion, su tamano ya es manejable (~840K aristas) y garantiza
    cobertura perfecta de los 6 tipos de nodo (las 7 tareas primarias
    ya tocan los 6 tipos entre todas).
  - GEN-cex-GEN (coexpressdb, 105.9M de las 110.3M aristas del grafo de
    entrenamiento, 96%): fraccion baja, 1.5%.
  - Resto de relaciones no primarias (~26 relaciones, ~3.19M aristas):
    fraccion mas alta, 25%, para preservar su señal.

Se aplica ademas la misma regla de conectividad de dos pasos que en la
Seccion 3: (1) toda arista con un extremo de grado global 1 (dentro del
grafo de entrenamiento, no del subgrafo) se incluye siempre; (2) tras
el muestreo, cualquier nodo que quede en grado 0 en el subgrafo se
rescata anadiendo de vuelta una de sus aristas del grafo de
entrenamiento completo.

Uso:
    python3 build_common_subgraph_v4.py \
        --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
        --out data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
        --report data/processed/evaluation/graph_v4/common_subgraph_report.json
"""
import argparse
import json

import numpy as np
import pandas as pd

PRIMARY_RELATIONS = {
    "GEN-ass-DIS", "GEN-ass-TIS", "GEN-ass-PWY", "CPD-trt-DIS",
    "CPD-int-GEN", "GEN-ppi-GEN", "CLL-mut-GEN",
}
COEXPRESSDB_RELATION = "GEN-cex-GEN"
FRAC_PRIMARY = 1.0
FRAC_COEXPRESSDB = 0.015
FRAC_OTHER = 0.25
SEED = 42


def log(msg):
    print(f"[subgraph] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    log("Cargando grafo de entrenamiento completo...")
    df = pd.read_csv(args.train_graph, sep="\t", dtype=str)
    df["node1_type"] = df["node1_type"].astype(str)
    df["node2_type"] = df["node2_type"].astype(str)
    log(f"{len(df):,} aristas, {df['relation'].nunique()} relaciones canonicas")

    log("Calculando grado global (sobre el grafo de entrenamiento completo)...")
    deg = pd.concat([df["node1"], df["node2"]]).value_counts()
    full_node_set = set(deg.index)
    log(f"{len(full_node_set):,} nodos unicos, grado minimo={deg.min()}, "
        f"nodos con grado 1: {(deg == 1).sum():,}")

    protected_mask = (deg.reindex(df["node1"]).to_numpy() == 1) | \
                      (deg.reindex(df["node2"]).to_numpy() == 1)
    n_protected = int(protected_mask.sum())
    log(f"Aristas protegidas (grado global 1 en algun extremo): {n_protected:,}")

    sampled_parts = []
    per_relation_report = []
    rng = np.random.default_rng(SEED)

    for relation, group in df.groupby("relation", observed=True):
        if relation in PRIMARY_RELATIONS:
            frac = FRAC_PRIMARY
            categoria = "primaria"
        elif relation == COEXPRESSDB_RELATION:
            frac = FRAC_COEXPRESSDB
            categoria = "coexpressdb"
        else:
            frac = FRAC_OTHER
            categoria = "otra"

        n_total_rel = len(group)
        if frac >= 1.0:
            sampled = group
        else:
            n_sample = max(1, int(round(n_total_rel * frac)))
            sampled = group.sample(n=n_sample, random_state=SEED)
        sampled_parts.append(sampled)
        per_relation_report.append({
            "relation": relation, "categoria": categoria,
            "n_total_train": int(n_total_rel), "fraccion_objetivo": frac,
            "n_muestreadas": int(len(sampled)),
        })

    sampled_df = pd.concat(sampled_parts, ignore_index=False)
    log(f"Aristas muestreadas por relacion (antes de proteccion/rescate): {len(sampled_df):,}")

    protected_df = df[protected_mask]
    combined = pd.concat([sampled_df, protected_df], ignore_index=False)
    combined = combined[~combined.index.duplicated(keep="first")]
    log(f"Aristas tras union con protegidas (dedup por indice): {len(combined):,}")

    subgraph_node_set = set(combined["node1"]) | set(combined["node2"])
    stranded = full_node_set - subgraph_node_set
    log(f"Nodos aislados tras muestreo (grado 0 en el subgrafo): {len(stranded):,}")

    n_rescued = 0
    if stranded:
        log("Rescatando nodos aislados (recuperando una arista del grafo de entrenamiento completo)...")
        # Prefiltro vectorizado: nos restringimos a las aristas que tocan
        # algun nodo aislado, en vez de iterar las 110.3M aristas completas.
        touch_mask = df["node1"].isin(stranded) | df["node2"].isin(stranded)
        candidates = df[touch_mask]
        log(f"  Aristas candidatas a rescate (tocan algun nodo aislado): {len(candidates):,}")

        node1_arr = candidates["node1"].to_numpy()
        node2_arr = candidates["node2"].to_numpy()
        cand_index = candidates.index.to_numpy()

        stranded_remaining = set(stranded)
        rescue_rows = []
        for i in range(len(candidates)):
            if not stranded_remaining:
                break
            a, b = node1_arr[i], node2_arr[i]
            need_a = a in stranded_remaining
            need_b = b in stranded_remaining
            if need_a or need_b:
                rescue_rows.append(cand_index[i])
                if need_a:
                    stranded_remaining.discard(a)
                if need_b:
                    stranded_remaining.discard(b)

        rescue_df = df.loc[df.index.isin(rescue_rows)]
        n_rescued = len(rescue_df)
        combined = pd.concat([combined, rescue_df], ignore_index=False)
        combined = combined[~combined.index.duplicated(keep="first")]
        log(f"Aristas rescatadas: {n_rescued}. Nodos que no se pudieron rescatar: {len(stranded_remaining)}")
        subgraph_node_set = set(combined["node1"]) | set(combined["node2"])

    combined = combined.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    combined.to_csv(args.out, sep="\t", index=False)
    log(f"Subgrafo guardado -> {args.out} ({len(combined):,} aristas, {len(subgraph_node_set):,} nodos)")

    # Cobertura por tipo de nodo, respecto a graph_v4 (Seccion 1 / MANIFEST_graph_v4.md)
    FULL_GRAPH_V4_NODE_COUNTS = {"GEN": 26330, "DIS": 23697, "CPD": 19629,
                                  "PWY": 2820, "CLL": 1953, "TIS": 699}
    node_type_map = {}
    for _, row in combined[["node1", "node1_type"]].drop_duplicates().iterrows():
        node_type_map[row["node1"]] = row["node1_type"]
    for _, row in combined[["node2", "node2_type"]].drop_duplicates().iterrows():
        node_type_map.setdefault(row["node2"], row["node2_type"])

    coverage_by_type = {}
    for t, total_v4 in FULL_GRAPH_V4_NODE_COUNTS.items():
        n_in_subgraph = sum(1 for v in node_type_map.values() if v == t)
        coverage_by_type[t] = {
            "n_graph_v4": total_v4, "n_subgrafo": n_in_subgraph,
            "cobertura_pct": round(100 * n_in_subgraph / total_v4, 2),
        }
        log(f"  Tipo {t}: {n_in_subgraph:,}/{total_v4:,} ({coverage_by_type[t]['cobertura_pct']}%)")

    report = {
        "seed": SEED,
        "metodo_muestreo": "estratificado por relacion canonica, con proteccion de "
                            "aristas de grado global 1 y rescate de nodos aislados residuales",
        "fracciones": {
            "relaciones_primarias": FRAC_PRIMARY,
            "GEN-cex-GEN_coexpressdb": FRAC_COEXPRESSDB,
            "resto_relaciones": FRAC_OTHER,
        },
        "n_aristas_grafo_entrenamiento": int(len(df)),
        "n_aristas_subgrafo": int(len(combined)),
        "fraccion_aristas_global": round(len(combined) / len(df), 5),
        "n_nodos_grafo_v4": len(full_node_set),
        "n_nodos_subgrafo": len(subgraph_node_set),
        "n_aristas_protegidas_grado1": n_protected,
        "n_nodos_rescatados": n_rescued,
        "cobertura_por_tipo_nodo": coverage_by_type,
        "detalle_por_relacion": per_relation_report,
    }
    with open(args.report, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    log(f"Reporte -> {args.report}")


if __name__ == "__main__":
    main()
