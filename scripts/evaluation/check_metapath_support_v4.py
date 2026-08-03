#!/usr/bin/env python3
"""
Seccion 5 (TFM graph_v4): verificacion de soporte de los 6 metapaths
propuestos para HIN2Vec, exigida explicitamente por el protocolo antes
de comprometerse a usarlos ("no incluyas un metapath por inercia").

Para cada metapath X-Y-X, se usa la relacion base X-rel-Y (o Y-rel-X)
en el grafo COMPLETO graph_v4 (no el split de train), se agrupan las
aristas por el nodo "hub" (el tipo Y del medio) y se mide:
  - cuantos nodos hub distintos hay
  - la distribucion de grado (cuantos X distintos toca cada hub)
  - cuantas instancias de metapath (pares X1-X2 no ordenados) implica
    cada hub -- C(grado,2) -- y que fraccion del total aporta el hub
    mas conectado (para detectar degeneracion tipo "todo conectado con
    todo" por un solo hub gigante, que no aporta señal especifica).

Metapaths evaluados (definicion, relacion base, tipo hub):
  GEN-DIS-GEN  <- GEN-ass-DIS   (hub=DIS)
  GEN-TIS-GEN  <- GEN-ass-TIS   (hub=TIS)
  GEN-PWY-GEN  <- GEN-ass-PWY   (hub=PWY)
  GEN-CPD-GEN  <- CPD-int-GEN   (hub=CPD, node1=CPD/node2=GEN)
  GEN-CLL-GEN  <- CLL-mut-GEN   (hub=CLL, node1=CLL/node2=GEN)
  CPD-DIS-CPD  <- CPD-trt-DIS   (hub=DIS, node1=CPD/node2=DIS)

Uso:
    python3 check_metapath_support_v4.py \
        --graph data/processed/integrated/graph_v4/graph_edges.tsv \
        --out data/processed/evaluation/graph_v4/metapath_support_v4.tsv
"""
import argparse
from math import comb

import pandas as pd

METAPATHS = [
    # (nombre_metapath, relacion_base, columna_hub_en_la_relacion)
    ("GEN-DIS-GEN", "GEN-ass-DIS", "node2"),  # hub = node2 (DIS)
    ("GEN-TIS-GEN", "GEN-ass-TIS", "node2"),  # hub = node2 (TIS)
    ("GEN-PWY-GEN", "GEN-ass-PWY", "node2"),  # hub = node2 (PWY)
    ("GEN-CPD-GEN", "CPD-int-GEN", "node1"),  # hub = node1 (CPD)
    ("GEN-CLL-GEN", "CLL-mut-GEN", "node1"),  # hub = node1 (CLL)
    ("CPD-DIS-CPD", "CPD-trt-DIS", "node2"),  # hub = node2 (DIS)
]


def log(msg):
    print(f"[metapath_support] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    relations_needed = {rel for _, rel, _ in METAPATHS}
    log(f"Relaciones base necesarias: {sorted(relations_needed)}")
    log(f"Leyendo {args.graph} en chunks, filtrando solo esas relaciones...")

    chunks_by_rel = {rel: [] for rel in relations_needed}
    n_rows_total = 0
    for chunk in pd.read_csv(args.graph, sep="\t", dtype=str,
                              usecols=["node1", "node2", "relation"],
                              chunksize=2_000_000):
        n_rows_total += len(chunk)
        sub = chunk[chunk["relation"].isin(relations_needed)]
        for rel, df in sub.groupby("relation"):
            chunks_by_rel[rel].append(df)

    log(f"{n_rows_total:,} aristas totales leidas de graph_v4")

    rows_out = []
    for name, rel, hub_col in METAPATHS:
        parts = chunks_by_rel.get(rel, [])
        if not parts:
            log(f"  {name}: relacion base '{rel}' NO encontrada en graph_v4 -- 0 soporte")
            rows_out.append({"metapath": name, "relacion_base": rel, "n_aristas_base": 0,
                              "n_hubs": 0, "grado_min": None, "grado_mediana": None,
                              "grado_max": None, "n_instancias_metapath": 0,
                              "frac_top_hub": None})
            continue
        df = pd.concat(parts, ignore_index=True)
        n_aristas_base = len(df)
        deg = df.groupby(hub_col).size()
        n_hubs = len(deg)
        n_instancias = int(sum(comb(int(d), 2) for d in deg))
        top_hub_instancias = comb(int(deg.max()), 2) if len(deg) else 0
        frac_top_hub = (top_hub_instancias / n_instancias) if n_instancias > 0 else 0.0

        log(f"  {name} (base={rel}, hub={hub_col}): {n_aristas_base:,} aristas, "
            f"{n_hubs:,} hubs distintos, grado[min/mediana/max]="
            f"{int(deg.min())}/{int(deg.median())}/{int(deg.max())}, "
            f"instancias_metapath={n_instancias:,}, frac_top_hub={frac_top_hub:.4f}")

        rows_out.append({
            "metapath": name, "relacion_base": rel, "n_aristas_base": n_aristas_base,
            "n_hubs": n_hubs, "grado_min": int(deg.min()), "grado_mediana": float(deg.median()),
            "grado_max": int(deg.max()), "n_instancias_metapath": n_instancias,
            "frac_top_hub": round(frac_top_hub, 6),
        })

    out_df = pd.DataFrame(rows_out)
    out_df.to_csv(args.out, sep="\t", index=False)
    log(f"Resultado -> {args.out}")


if __name__ == "__main__":
    main()
