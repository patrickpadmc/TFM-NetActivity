#!/usr/bin/env python
"""
Convierte el archivo .emb generado por pecanpy (formato texto gensim
word2vec: primera linea "n_nodes dim", luego "node_id v1 v2 ... vd")
a un .npy + node_order.tsv, siguiendo la misma convencion de salida
que spectral_embedding_v3.py.

Parseo manual linea por linea en vez de pd.read_csv: algunos node IDs
del grafo contienen espacios (rompe el conteo de columnas de un parser
CSV estandar), asi que se toman los ultimos `dim` tokens como el vector
y se reconstruye el node ID con todo lo anterior.

Uso:
    python convert_embeddings_v3.py --emb <path al .emb de pecanpy> --out-dir <path> --dim 256
"""
import argparse
import os
import time
import numpy as np
import pandas as pd


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--emb", required=True, help="Path al .emb generado por pecanpy")
    parser.add_argument("--out-dir", required=True, help="Directorio de salida")
    parser.add_argument("--dim", type=int, default=256, help="Dimension del embedding")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    dim = args.dim

    log(f"Leyendo {args.emb}")
    node_ids = []
    vectors = []
    with open(args.emb, "r") as f:
        header = f.readline()
        n_nodes_expected, dim_expected = header.split()
        n_nodes_expected = int(n_nodes_expected)
        dim_expected = int(dim_expected)
        if dim_expected != dim:
            raise ValueError(f"Dimension en header ({dim_expected}) != --dim ({dim})")

        for line_num, line in enumerate(f, start=2):
            parts = line.rstrip("\n").split(" ")
            if len(parts) < dim + 1:
                raise ValueError(f"Linea {line_num}: se esperaban al menos {dim + 1} campos, "
                                  f"se encontraron {len(parts)}")
            vec_tokens = parts[-dim:]
            id_tokens = parts[:-dim]
            node_id = " ".join(id_tokens)
            node_ids.append(node_id)
            vectors.append([float(x) for x in vec_tokens])

    embeddings = np.asarray(vectors, dtype=np.float32)
    log(f"Cargados {embeddings.shape[0]} nodos, dim={embeddings.shape[1]} "
        f"(header decia {n_nodes_expected} nodos)")

    if embeddings.shape[0] != n_nodes_expected:
        log(f"ADVERTENCIA: el header decia {n_nodes_expected} nodos pero se parsearon "
            f"{embeddings.shape[0]}")

    node_order_path = os.path.join(args.out_dir, "node_order_node2vec_v3.tsv")
    pd.Series(node_ids).to_csv(node_order_path, sep="\t", index=False, header=False)
    log(f"Guardado orden de nodos en {node_order_path}")

    emb_path = os.path.join(args.out_dir, f"node2vec_embeddings_v3_d{dim}.npy")
    np.save(emb_path, embeddings)
    log(f"Embeddings guardados {embeddings.shape} en {emb_path}")

    log("Listo")


if __name__ == "__main__":
    main()
