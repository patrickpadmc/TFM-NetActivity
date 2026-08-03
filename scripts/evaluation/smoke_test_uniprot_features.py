#!/usr/bin/env python3
"""
Seccion 10: smoke test rapido (CPU, datos sinteticos minusculos) de la
modificacion de gae_common.py / gat_common.py para las condiciones
uniprot_fixed / uniprot_finetuned, ANTES de lanzar las 4 corridas reales
en GPU. Corre en segundos. No usa datos reales del grafo.

Verifica:
  1. build_graph + load_uniprot_features + NodeFeatureModule encajan sin
     errores de forma/dtype para GAE y GAT-AE.
  2. 'structural' (uniprot_tsv=None) sigue funcionando exactamente igual
     que antes (no se rompe nada de lo ya validado en Secciones 4-9).
  3. 'fixed': uniprot_raw.requires_grad es False y sus valores NO cambian
     tras un paso de entrenamiento.
  4. 'finetuned': uniprot_raw.requires_grad es True y sus valores SI
     cambian tras un paso de entrenamiento.
  5. Los nodos GEN sin cobertura (no listados en el tsv de prueba) y los
     nodos no-GEN siguen recibiendo gradiente en su fila de
     nn.Embedding de respaldo (no se quedan congelados por accidente).

Uso (en un nodo de computo, via srun/sbatch -- nunca en el login node):
    cd ${REPO}/scripts/evaluation
    python3 smoke_test_uniprot_features.py
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gae_common import build_graph, load_uniprot_features, GAEEncoder, train_gae
from gat_common import GATEncoder, train_gat


def _make_tiny_graph(tmpdir, n_nodes=24, n_edges=80, seed=0):
    rng = np.random.default_rng(seed)
    # la mitad de los nodos se nombran como si fueran GEN (ENSG...), la otra mitad no
    nodes = [f"ENSG_{i}" for i in range(n_nodes // 2)] + [f"DIS_{i}" for i in range(n_nodes // 2)]
    edges = set()
    while len(edges) < n_edges:
        u, v = rng.choice(nodes, size=2, replace=False)
        if u != v:
            edges.add((u, v))
    df = pd.DataFrame(list(edges), columns=["node1", "node2"])
    path = os.path.join(tmpdir, "tiny_graph.tsv")
    df.to_csv(path, sep="\t", index=False)
    return path, nodes


def _make_tiny_uniprot(tmpdir, covered_gen_nodes, ext_dim=8, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for n in covered_gen_nodes:
        vec = rng.normal(size=ext_dim).astype(np.float32)
        rows.append([n, "GEN", "UniProt_ProtT5", ext_dim, *vec])
    cols = ["node_id", "node_type", "embedding_source", "dimension"] + [f"e_{i+1:04d}" for i in range(ext_dim)]
    df = pd.DataFrame(rows, columns=cols)
    path = os.path.join(tmpdir, "tiny_uniprot.tsv")
    df.to_csv(path, sep="\t", index=False)
    return path


def check_load_uniprot_features(tmpdir):
    print("== 1. load_uniprot_features + alineacion con idx ==")
    graph_path, nodes = _make_tiny_graph(tmpdir)
    A, pos_pairs, nodes2, idx = build_graph(graph_path)
    assert set(nodes2) <= set(nodes), "build_graph produjo nodos inesperados"
    gen_nodes_in_graph = [n for n in nodes2 if n.startswith("ENSG_")]
    covered = gen_nodes_in_graph[: max(1, len(gen_nodes_in_graph) // 2)]  # cubrir solo la mitad
    uniprot_path = _make_tiny_uniprot(tmpdir, covered)
    uniprot_node_indices, uniprot_vectors = load_uniprot_features(uniprot_path, idx)
    assert len(uniprot_node_indices) == len(covered), \
        f"esperaba {len(covered)} nodos cubiertos, obtuve {len(uniprot_node_indices)}"
    assert uniprot_vectors.shape == (len(covered), 8)
    recovered_ids = {nodes2[i] for i in uniprot_node_indices}
    assert recovered_ids == set(covered), "los indices no corresponden a los node_id esperados"
    print(f"   OK: {len(covered)} nodos GEN cubiertos de {len(nodes2)} totales, alineados correctamente.")
    return graph_path, uniprot_path, idx, nodes2, covered


def check_encoder_structural(EncoderClass, forward_args_fn, n_nodes, label):
    print(f"== 2. {label}: condicion 'structural' (sin UniProt) sigue funcionando ==")
    model = EncoderClass(n_nodes, in_dim=16, **({"hidden_dim": 16, "out_dim": 16} if EncoderClass is GAEEncoder
                                                  else {"hidden_dim": 8, "heads1": 2, "out_dim": 16}))
    assert model.feature_module.has_uniprot is False
    out = model(*forward_args_fn(model))
    assert out.shape == (n_nodes, 16)
    assert torch.isfinite(out).all()
    print(f"   OK: forward sin errores, salida ({n_nodes}, 16), sin NaN/Inf.")


def check_uniprot_mode(EncoderClass, forward_args_fn, n_nodes, uniprot_vectors, uniprot_node_indices,
                        mode, label):
    print(f"== 3. {label}: condicion 'uniprot_{mode}' ==")
    kwargs = {"hidden_dim": 16, "out_dim": 16} if EncoderClass is GAEEncoder else \
             {"hidden_dim": 8, "heads1": 2, "out_dim": 16}
    model = EncoderClass(n_nodes, in_dim=16, uniprot_vectors=uniprot_vectors,
                          uniprot_node_indices=uniprot_node_indices, uniprot_mode=mode, **kwargs)
    fm = model.feature_module
    assert fm.has_uniprot is True
    if mode == "fixed":
        assert fm.uniprot_raw.requires_grad is False, "fixed deberia congelar uniprot_raw"
    else:
        assert fm.uniprot_raw.requires_grad is True, "finetuned deberia permitir gradiente en uniprot_raw"

    initial_snapshot = fm.initial_projection_snapshot()
    assert initial_snapshot.shape == (len(uniprot_node_indices), 16)
    raw_before = fm.uniprot_raw.detach().clone()
    backup_emb_before = fm.node_features.weight.detach().clone()

    out = model(*forward_args_fn(model))
    assert out.shape == (n_nodes, 16)
    assert torch.isfinite(out).all()

    loss = out.sum()
    loss.backward()

    opt = torch.optim.SGD(model.parameters(), lr=0.1)
    opt.step()

    raw_after = fm.uniprot_raw.detach().clone()
    backup_emb_after = fm.node_features.weight.detach().clone()

    if mode == "fixed":
        assert torch.allclose(raw_before, raw_after), \
            "fixed: uniprot_raw NO deberia cambiar tras un paso de entrenamiento, pero cambio"
    else:
        assert not torch.allclose(raw_before, raw_after), \
            "finetuned: uniprot_raw DEBERIA cambiar tras un paso de entrenamiento, pero no cambio"

    # los nodos de respaldo (no cubiertos por UniProt) deben seguir recibiendo gradiente
    uncovered_mask = torch.ones(n_nodes, dtype=torch.bool)
    uncovered_mask[uniprot_node_indices] = False
    assert not torch.allclose(backup_emb_before[uncovered_mask], backup_emb_after[uncovered_mask]), \
        "los nodos de respaldo (no cubiertos) deberian seguir aprendiendo via nn.Embedding"

    print(f"   OK: requires_grad correcto, forward sin errores, "
          f"uniprot_raw {'NO cambia (correcto)' if mode == 'fixed' else 'SI cambia (correcto)'}, "
          f"nodos de respaldo siguen aprendiendo.")


def check_train_loop(train_fn, graph_path, uniprot_path, mode, label):
    print(f"== 4. {label}: train_{'gae' if train_fn is train_gae else 'gat'} end-to-end, uniprot_{mode} ==")
    kwargs = dict(dim=16, max_epochs=2, patience=100, check_every=1,
                  n_val_internal=2, seed=42, device=torch.device("cpu"),
                  uniprot_tsv=uniprot_path, uniprot_mode=mode)
    if train_fn is train_gae:
        kwargs.update(hidden_dim=16, batch_edges=20)
    else:
        kwargs.update(hidden_dim=8, heads1=2, batch_edges=20, use_amp=False)
    embeddings, timings = train_fn(graph_path, **kwargs)
    assert timings["uniprot_mode"] == mode
    assert timings["n_uniprot_covered"] > 0
    assert timings["initial_uniprot_embeddings"] is not None
    n_gen_check = list(timings["initial_uniprot_embeddings"].keys())[0]
    assert n_gen_check in embeddings
    print(f"   OK: {timings['n_uniprot_covered']} nodos cubiertos, "
          f"initial_uniprot_embeddings presente ({len(timings['initial_uniprot_embeddings'])} vectores), "
          f"embeddings finales generados para {len(embeddings)} nodos.")


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        graph_path, uniprot_path, idx, nodes2, covered = check_load_uniprot_features(tmpdir)
        n_nodes = len(nodes2)
        uniprot_node_indices, uniprot_vectors = load_uniprot_features(uniprot_path, idx)

        # Para GAE necesitamos A_norm real; reconstruimos con normalized_adj_torch.
        from gae_common import normalized_adj_torch
        A, _, _, _ = build_graph(graph_path)
        A_norm = normalized_adj_torch(A, torch.device("cpu"))
        gae_fwd = lambda m: (A_norm,)

        from gat_common import build_edge_index
        edge_index = build_edge_index(A, torch.device("cpu"))
        gat_fwd = lambda m: (edge_index, n_nodes)

        check_encoder_structural(GAEEncoder, gae_fwd, n_nodes, "GAE")
        check_encoder_structural(GATEncoder, gat_fwd, n_nodes, "GAT-AE")

        for mode in ("fixed", "finetuned"):
            check_uniprot_mode(GAEEncoder, gae_fwd, n_nodes, uniprot_vectors, uniprot_node_indices, mode, "GAE")
            check_uniprot_mode(GATEncoder, gat_fwd, n_nodes, uniprot_vectors, uniprot_node_indices, mode, "GAT-AE")

        for mode in ("fixed", "finetuned"):
            check_train_loop(train_gae, graph_path, uniprot_path, mode, "GAE")
            check_train_loop(train_gat, graph_path, uniprot_path, mode, "GAT-AE")

        print("\nTODOS LOS CHEQUEOS PASARON.")


if __name__ == "__main__":
    main()
