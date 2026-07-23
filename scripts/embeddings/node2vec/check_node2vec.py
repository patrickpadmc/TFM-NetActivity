import numpy as np
import pandas as pd

out_dir = "data/processed/embeddings/node2vec_v3"

emb = np.load(f"{out_dir}/node2vec_embeddings_v3_d256.npy")
node_order = pd.read_csv(f"{out_dir}/node_order_node2vec_v3.tsv", sep="\t", header=None)[0]

print("shape embeddings:", emb.shape)
print("shape node_order:", node_order.shape)
print("nodos duplicados en node_order:", node_order.duplicated().sum())
print("NaNs:", np.isnan(emb).sum(), "| Infs:", np.isinf(emb).sum())

norms = np.linalg.norm(emb, axis=1)
print()
print("norma: min={:.4f} max={:.4f} media={:.4f} std={:.4f}".format(
    norms.min(), norms.max(), norms.mean(), norms.std()))

print()
print("media por dimension (primeras 5):", emb.mean(axis=0)[:5])
print("std por dimension (primeras 5):", emb.std(axis=0)[:5])
