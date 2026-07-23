#!/usr/bin/env python
"""
Lee graph_v3/graph_edges.tsv una sola vez y cachea a parquet.

IMPORTANTE: la columna `relation` del tsv ya viene en formato
"node1_type-codigo-node2_type" (ej. "GEN-ass-DIS", "CPD-upe-GEN"), NO
es un codigo corto ("ass") como sugeria el nombre de las carpetas en
los logs de carga del grafo. Por eso edge_class = relation, sin
reconstruir nada -- concatenar node1_type + relation + node2_type
duplicaba los tipos de nodo (bug encontrado en produccion: daba
"GEN-GEN-ass-DIS-DIS" en vez de "GEN-ass-DIS").

Uso:
    python prepare_edge_class_v3.py --edges <path a graph_edges.tsv> --out <path al parquet cacheado>
"""
import argparse
import time
import pyarrow.csv as pv
import pyarrow.parquet as pq


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", required=True, help="Path a graph_edges.tsv de graph_v3")
    parser.add_argument("--out", required=True, help="Path del parquet cacheado de salida")
    args = parser.parse_args()

    log(f"Leyendo {args.edges}")
    parse_opts = pv.ParseOptions(delimiter="\t")
    convert_opts = pv.ConvertOptions(
        include_columns=["node1", "node2", "node1_type", "node2_type", "relation"]
    )
    table = pv.read_csv(args.edges, parse_options=parse_opts, convert_options=convert_opts)
    log(f"{table.num_rows} filas cargadas")

    log("edge_class = relation (ya viene en formato node1_type-codigo-node2_type)")
    table = table.append_column("edge_class", table.column("relation"))
    log(f"Ejemplos de edge_class: {table.column('edge_class')[:3].to_pylist()}")

    log(f"Escribiendo {args.out}")
    pq.write_table(table, args.out)
    log("Listo")


if __name__ == "__main__":
    main()
