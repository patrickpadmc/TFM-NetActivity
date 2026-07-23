#!/usr/bin/env python
"""
Entrena HIN2Vec sobre una variante de graph_v3 (edges ya filtrados por
build_hin_input_v3.py) y guarda embeddings de nodo, embeddings de
metapath (relacion-secuencia) y el modelo entrenado.

Hiperparametros por defecto: embed_size=256, walk_length=100, walk=10,
window=4, neg=5, epochs=10, batch_size=128 (confirmados).

Uso:
    python train_hin2vec_v3.py --edges-parquet <path> --out-dir <path> --run-name auto_cex
"""
import argparse
import os
import pickle
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor"))
from walker import load_a_HIN_from_pandas
from model import NSTrainSet, HIN2vec, train


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges-parquet", required=True, help="Salida de build_hin_input_v3.py")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--run-name", required=True,
                         help="ej. auto_cex, auto_nocex, metapath_cex, metapath_nocex")
    parser.add_argument("--window", type=int, default=4)
    parser.add_argument("--walk", type=int, default=10, help="walks por nodo de partida")
    parser.add_argument("--walk-length", type=int, default=100)
    parser.add_argument("--embed-size", type=int, default=256)
    parser.add_argument("--neg", type=int, default=5)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--max-samples", type=int, default=30_000_000,
                         help="Cap de muestras positivas tras hin.sample() (submuestreo con seed fija)")
    parser.add_argument("--sample-seed", type=int, default=42)
    parser.add_argument("--num-workers", type=int, default=8)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    log(f"device = {device}")

    log(f"Leyendo {args.edges_parquet}")
    edges = pd.read_parquet(args.edges_parquet)
    log(f"{len(edges)} aristas (una direccion) cargadas")

    log("Construyendo HIN (grafo interno de walker.py)")
    hin = load_a_HIN_from_pandas([edges])
    hin.window = args.window
    log(f"Grafo: {hin.node_size} nodos")

    log(f"Generando walks (length={args.walk_length}, walks/nodo={args.walk}) "
        f"y muestreando pares (start, end, path_id)")
    t0 = time.time()
    samples = hin.sample(args.walk_length, args.walk,
                          max_samples=args.max_samples, seed=args.sample_seed)
    log(f"{len(samples)} muestras conservadas en {time.time() - t0:.1f}s, "
        f"path_size (relaciones distintas encontradas)={hin.path_size}")

    dataset = NSTrainSet(samples, hin.node_size, neg=args.neg)
    hin2vec = HIN2vec(hin.node_size, hin.path_size, args.embed_size, sigmoid_reg=True).to(device)

    data_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True,
                              num_workers=args.num_workers, pin_memory=True)
    optimizer = optim.AdamW(hin2vec.parameters())
    loss_function = nn.BCELoss()

    log(f"Entrenando {args.epochs} epochs")
    for epoch in range(args.epochs):
        train(200, hin2vec, device, data_loader, optimizer, loss_function, epoch)

    log("Entrenamiento completo, guardando outputs")

    node_embeds = hin2vec.start_embeds.weight.data.cpu().numpy()
    node_names = [hin.id2node[i] for i in range(hin.node_size)]

    emb_path = os.path.join(args.out_dir, f"hin2vec_embeddings_v3_{args.run_name}_d{args.embed_size}.npy")
    order_path = os.path.join(args.out_dir, f"node_order_hin2vec_v3_{args.run_name}.tsv")
    np.save(emb_path, node_embeds)
    pd.Series(node_names).to_csv(order_path, sep="\t", index=False, header=False)
    log(f"Node embeddings guardados {node_embeds.shape} en {emb_path}")

    path_embeds = hin2vec.path_embeds.weight.data.cpu().numpy()
    path_map = {pid: "|".join(rel_seq) for rel_seq, pid in hin._path2id.items()}
    path_emb_path = os.path.join(args.out_dir, f"hin2vec_path_embeddings_v3_{args.run_name}.npy")
    path_map_path = os.path.join(args.out_dir, f"path_id_map_v3_{args.run_name}.pkl")
    np.save(path_emb_path, path_embeds)
    with open(path_map_path, "wb") as f:
        pickle.dump(path_map, f)
    log(f"Path embeddings guardados {path_embeds.shape} en {path_emb_path}")
    log(f"Mapa path_id -> secuencia de relaciones guardado en {path_map_path}")

    model_path = os.path.join(args.out_dir, f"hin2vec_model_v3_{args.run_name}.pt")
    torch.save(hin2vec.state_dict(), model_path)
    log(f"Modelo guardado en {model_path}")

    log("Listo")


if __name__ == "__main__":
    main()
