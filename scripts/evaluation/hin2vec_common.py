#!/usr/bin/env python3
"""
Seccion 5 (TFM graph_v4): HIN2Vec (Fu, Lee, Lei 2017), implementado a
mano en PyTorch puro, siguiendo el patron ya usado en GAE/GAT-AE
(sin paquete disponible, decision confirmada por el usuario).

DISEÑO (confirmado por el usuario antes de implementar):
  - Caminatas restringidas a metapath: para cada uno de los 6 metapaths
    confirmados (todos de la forma X-Y-X, 2 saltos, simetricos), se
    camina alternando X -> hub(Y) -> X -> hub(Y) -> ... eligiendo
    vecinos uniformemente al azar en cada salto (sin sesgo p,q tipo
    Node2Vec -- la restriccion de tipo/metapath ya es la innovacion
    central de HIN2Vec, no hace falta sesgo adicional).
  - De cada caminata se extraen SOLO pares de nodos X adyacentes
    (posiciones consecutivas en la secuencia de nodos X de la
    caminata, ej. x0-x1, x1-x2, ...) como instancias positivas del
    metapath -- NO se usa ventana tipo skip-gram mas amplia, porque
    eso mezclaria metapaths compuestos mas largos (X-Y-X-Y-X) con el
    metapath de 2 saltos pedido explicitamente por el protocolo.
  - Modelo: embedding de nodo compartido entre todos los tipos (dim
    128, igual que GAE/GAT-AE, sin atributos biologicos) + embedding
    de metapath/relacion (dim 128, uno por cada uno de los 6
    metapaths). Dado (u, v, r): compuerta = sigmoid(emb_relacion(r))
    (regularizacion tipo "auto-atencion binaria" del paper original,
    fuerza a los pesos de la relacion a comportarse como selectores
    de features ~binarios); h = emb_nodo(u) * emb_nodo(v) * compuerta;
    logit = sum(h). Entrenado como clasificacion binaria (BCE) contra
    negativos.
  - Negative sampling: se corrompe v (o u) por un nodo X aleatorio,
    SIN verificar exhaustivamente que no sea un verdadero vecino de 2
    saltos (verificarlo exigiria materializar la adyacencia X-X
    completa, exactamente lo que se evito en el chequeo de soporte de
    metapaths por el riesgo de explosion combinatoria en hubs grandes).
    Es el mismo tipo de aproximacion que usa word2vec/Node2Vec con
    muestreo por unigrama -- una pequeña tasa de falsos negativos es
    aceptada y documentada, no corregida.
  - GEN-TIS-GEN usa el filtro de grado confirmado (excluir tejidos con
    grado > 1147, percentil 90 de su propia distribucion) SOLO para
    elegir que nodos TIS actuan como hub en las caminatas -- todos los
    nodos TIS igual reciben una fila de embedding entrenable (para no
    perjudicar la evaluacion de GEN-ass-TIS en la Seccion 4/5).

CARACTERISTICAS DE NODO: identico a GAE/GAT-AE -- embedding entrenable
por nodo, sin atributos biologicos.
"""
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from numba import njit, prange
from sklearn.metrics import average_precision_score

# (nombre_metapath, relacion_base, columna_X_en_la_relacion,
#  columna_hub_en_la_relacion, grado_maximo_hub_o_None)
METAPATH_DEF = [
    ("GEN-DIS-GEN", "GEN-ass-DIS", "node1", "node2", None),
    ("GEN-TIS-GEN", "GEN-ass-TIS", "node1", "node2", 1147),
    ("GEN-PWY-GEN", "GEN-ass-PWY", "node1", "node2", None),
    ("GEN-CPD-GEN", "CPD-int-GEN", "node2", "node1", None),
    ("GEN-CLL-GEN", "CLL-mut-GEN", "node2", "node1", None),
    ("CPD-DIS-CPD", "CPD-trt-DIS", "node1", "node2", None),
]
METAPATH_NAMES = [m[0] for m in METAPATH_DEF]


def log(msg):
    print(f"[hin2vec] {msg}", flush=True)


def _build_csr(pairs, n_left, n_right):
    """pairs: array (M,2) de indices [left, right]. Devuelve CSR left->right."""
    order = np.argsort(pairs[:, 0], kind="stable")
    pairs_sorted = pairs[order]
    indptr = np.zeros(n_left + 1, dtype=np.int64)
    counts = np.bincount(pairs_sorted[:, 0], minlength=n_left)
    indptr[1:] = np.cumsum(counts)
    indices = pairs_sorted[:, 1].astype(np.int64)
    return indptr, indices


def build_all_metapath_data(edges_path):
    """
    Carga el TSV de aristas (node1,node2,node1_type,node2_type,relation,...)
    y construye, para cada uno de los 6 metapaths: la adyacencia bipartita
    X<->hub (con el filtro de grado aplicado solo al lado hub cuando
    corresponda) usando un UNICO espacio de indices de nodo global
    (compartido entre todos los metapaths y todos los tipos de nodo).

    Devuelve:
      node_list: lista de node_id (str) ordenada, indice global -> id
      node_index: dict id -> indice global
      metapath_data: dict nombre -> dict con:
          indptr_x, indices_x (X -> hub, indices locales a X y a hub)
          indptr_hub, indices_hub (hub -> X)
          x_global_ids: array de indices GLOBALES para cada indice local X
          hub_global_ids: array de indices GLOBALES para cada indice local hub
          x_local_of_global: dict indice_global -> indice_local_X (para arrancar caminatas)
    """
    needed_rel = {rel for _, rel, _, _, _ in METAPATH_DEF}
    log(f"Cargando {edges_path}, filtrando relaciones: {sorted(needed_rel)}")
    chunks = []
    for chunk in pd.read_csv(edges_path, sep="\t", dtype=str,
                              usecols=["node1", "node2", "relation"],
                              chunksize=2_000_000):
        sub = chunk[chunk["relation"].isin(needed_rel)]
        if len(sub):
            chunks.append(sub)
    df_all = pd.concat(chunks, ignore_index=True)
    log(f"{len(df_all):,} aristas relevantes cargadas")

    all_nodes = pd.unique(pd.concat([df_all["node1"], df_all["node2"]], ignore_index=True))
    node_list = sorted(all_nodes.tolist())
    node_index = {nid: i for i, nid in enumerate(node_list)}
    log(f"{len(node_list):,} nodos distintos en el espacio global de embeddings")

    metapath_data = {}
    for name, rel, x_col, hub_col, max_deg in METAPATH_DEF:
        df = df_all[df_all["relation"] == rel]
        if len(df) == 0:
            log(f"  {name}: sin aristas para la relacion base '{rel}', se omite")
            continue

        x_ids = df[x_col].values
        hub_ids = df[hub_col].values

        # filtro de grado (solo afecta que nodos hub sirven para caminar,
        # no remueve filas del embedding global)
        if max_deg is not None:
            hub_deg = pd.Series(hub_ids).value_counts()
            hubs_ok = set(hub_deg[hub_deg <= max_deg].index)
            mask = np.array([h in hubs_ok for h in hub_ids])
            n_excluidos = len(set(hub_ids.tolist())) - len(hubs_ok)
            log(f"  {name}: filtro de grado <= {max_deg} excluye {n_excluidos} hubs "
                f"({mask.sum():,}/{len(mask):,} aristas retenidas)")
            x_ids, hub_ids = x_ids[mask], hub_ids[mask]

        x_unique = sorted(set(x_ids.tolist()))
        hub_unique = sorted(set(hub_ids.tolist()))
        x_local_index = {nid: i for i, nid in enumerate(x_unique)}
        hub_local_index = {nid: i for i, nid in enumerate(hub_unique)}

        pairs_xh = np.array([[x_local_index[x], hub_local_index[h]]
                              for x, h in zip(x_ids, hub_ids)], dtype=np.int64)
        pairs_hx = pairs_xh[:, [1, 0]]

        indptr_x, indices_x = _build_csr(pairs_xh, len(x_unique), len(hub_unique))
        indptr_hub, indices_hub = _build_csr(pairs_hx, len(hub_unique), len(x_unique))

        x_global_ids = np.array([node_index[nid] for nid in x_unique], dtype=np.int64)
        hub_global_ids = np.array([node_index[nid] for nid in hub_unique], dtype=np.int64)

        log(f"  {name}: {len(x_unique):,} nodos X, {len(hub_unique):,} nodos hub, "
            f"{len(pairs_xh):,} aristas base usadas para caminar")

        metapath_data[name] = {
            "indptr_x": indptr_x, "indices_x": indices_x,
            "indptr_hub": indptr_hub, "indices_hub": indices_hub,
            "x_global_ids": x_global_ids, "hub_global_ids": hub_global_ids,
            "n_x": len(x_unique), "n_hub": len(hub_unique),
        }

    return node_list, node_index, metapath_data


@njit(cache=True, parallel=True)
def _generate_bipartite_walks(indptr_x, indices_x, indptr_hub, indices_hub,
                                start_x_nodes, walk_length, seed):
    """Caminatas X -> hub -> X -> hub -> ... arrancando desde CADA nodo X
    (start_x_nodes = np.arange(n_x)). Devuelve indices LOCALES a X en cada
    posicion par de la caminata (0, 2, 4, ...); -1 si la caminata se corto
    antes.

    Nota: cur_x se inicializa desde start_x_nodes[i] (un array), NO
    directamente desde el indice de bucle i -- si se hace cur_x = i y
    luego se reasigna cur_x dentro del loop, el pase parfor de numba lo
    interpreta como "overwrite of parallel loop index" y falla la
    compilacion (visto empiricamente). Desacoplar cur_x de i via un array
    evita la confusion del analisis de numba, igual que en node2vec_common.py.
    """
    n_starts = len(start_x_nodes)
    n_x_positions = (walk_length + 1) // 2
    walks = np.full((n_starts, n_x_positions), -1, dtype=np.int64)
    for i in prange(n_starts):
        np.random.seed(seed + i)
        cur_x = start_x_nodes[i]
        walks[i, 0] = cur_x
        pos = 1
        for step in range(walk_length):
            xs, xe = indptr_x[cur_x], indptr_x[cur_x + 1]
            if xe <= xs:
                break
            hub = indices_x[xs + np.random.randint(0, xe - xs)]
            hs, he = indptr_hub[hub], indptr_hub[hub + 1]
            if he <= hs:
                break
            nxt_x = indices_hub[hs + np.random.randint(0, he - hs)]
            if step % 2 == 1:
                walks[i, pos] = nxt_x
                pos += 1
            cur_x = nxt_x
    return walks


def generate_metapath_pairs(mp_info, walk_length=6, seed=42):
    """Genera pares (u,v) GLOBALES a partir de caminatas restringidas al
    metapath, tomando solo posiciones X consecutivas de cada caminata."""
    start_x_nodes = np.arange(mp_info["n_x"], dtype=np.int64)
    walks_local = _generate_bipartite_walks(
        mp_info["indptr_x"], mp_info["indices_x"],
        mp_info["indptr_hub"], mp_info["indices_hub"],
        start_x_nodes, walk_length, seed,
    )
    u_list, v_list = [], []
    for col in range(walks_local.shape[1] - 1):
        u_loc = walks_local[:, col]
        v_loc = walks_local[:, col + 1]
        valid = (u_loc >= 0) & (v_loc >= 0)
        u_list.append(u_loc[valid])
        v_list.append(v_loc[valid])
    u_loc_all = np.concatenate(u_list) if u_list else np.zeros(0, dtype=np.int64)
    v_loc_all = np.concatenate(v_list) if v_list else np.zeros(0, dtype=np.int64)
    u_glob = mp_info["x_global_ids"][u_loc_all]
    v_glob = mp_info["x_global_ids"][v_loc_all]
    return u_glob, v_glob


class HIN2VecModel(nn.Module):
    def __init__(self, n_nodes, n_metapaths, dim=128):
        super().__init__()
        self.node_emb = nn.Embedding(n_nodes, dim)
        self.rel_emb = nn.Embedding(n_metapaths, dim)
        nn.init.xavier_uniform_(self.node_emb.weight)
        nn.init.xavier_uniform_(self.rel_emb.weight)

    def forward(self, u_idx, v_idx, r_idx):
        u = self.node_emb(u_idx)
        v = self.node_emb(v_idx)
        r_gate = torch.sigmoid(self.rel_emb(r_idx))
        h = u * v * r_gate
        return h.sum(dim=-1)


def train_hin2vec(edges_path, dim=128, walk_length=6, seed=42, lr=0.005,
                   weight_decay=0.0, max_epochs=200, patience=20, check_every=5,
                   batch_size=200_000, n_val_per_metapath=2_000, device=None,
                   log_prefix="hin2vec"):
    def _log(msg):
        print(f"[{log_prefix}] {msg}", flush=True)

    torch.manual_seed(seed)
    np.random.seed(seed)
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Device: {device}")

    t0 = time.time()
    node_list, node_index, metapath_data = build_all_metapath_data(edges_path)
    n_nodes = len(node_list)
    metapaths_presentes = [m for m in METAPATH_NAMES if m in metapath_data]
    mp_id = {m: i for i, m in enumerate(metapaths_presentes)}
    _log(f"Metapaths con datos: {metapaths_presentes}. t={time.time()-t0:.1f}s")

    # generar pares positivos una vez por metapath (caminatas), separar
    # una porcion pequeña para validacion interna (early stopping)
    train_u, train_v, train_r = [], [], []
    val_u, val_v, val_r = [], [], []
    rng = np.random.default_rng(seed)
    for name in metapaths_presentes:
        u_glob, v_glob = generate_metapath_pairs(metapath_data[name], walk_length=walk_length, seed=seed)
        n_pairs = len(u_glob)
        n_val = min(n_val_per_metapath, max(1, n_pairs // 10))
        perm = rng.permutation(n_pairs)
        val_idx, train_idx = perm[:n_val], perm[n_val:]
        train_u.append(u_glob[train_idx]); train_v.append(v_glob[train_idx])
        train_r.append(np.full(len(train_idx), mp_id[name], dtype=np.int64))
        val_u.append(u_glob[val_idx]); val_v.append(v_glob[val_idx])
        val_r.append(np.full(len(val_idx), mp_id[name], dtype=np.int64))
        _log(f"  {name}: {n_pairs:,} pares positivos (caminatas), {n_val:,} para validacion interna")

    train_u = np.concatenate(train_u); train_v = np.concatenate(train_v); train_r = np.concatenate(train_r)
    val_u = np.concatenate(val_u); val_v = np.concatenate(val_v); val_r = np.concatenate(val_r)
    n_train = len(train_u)
    _log(f"Total: {n_train:,} pares de entrenamiento, {len(val_u):,} de validacion interna. t={time.time()-t0:.1f}s")

    val_neg_v = rng.integers(0, n_nodes, size=len(val_u))
    val_u_t = torch.tensor(np.concatenate([val_u, val_u]), dtype=torch.long, device=device)
    val_v_t = torch.tensor(np.concatenate([val_v, val_neg_v]), dtype=torch.long, device=device)
    val_r_t = torch.tensor(np.concatenate([val_r, val_r]), dtype=torch.long, device=device)
    val_labels = np.concatenate([np.ones(len(val_u)), np.zeros(len(val_u))])

    model = HIN2VecModel(n_nodes, len(metapaths_presentes), dim=dim).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    loss_fn = nn.BCEWithLogitsLoss()

    best_val_auprc, best_state, epochs_no_improve = -1.0, None, 0
    history = []

    for epoch in range(1, max_epochs + 1):
        te0 = time.time()
        model.train()
        batch_idx = rng.choice(n_train, size=min(batch_size, n_train), replace=False)
        bu, bv, br = train_u[batch_idx], train_v[batch_idx], train_r[batch_idx]
        neg_v = rng.integers(0, n_nodes, size=len(bu))

        u_t = torch.tensor(np.concatenate([bu, bu]), dtype=torch.long, device=device)
        v_t = torch.tensor(np.concatenate([bv, neg_v]), dtype=torch.long, device=device)
        r_t = torch.tensor(np.concatenate([br, br]), dtype=torch.long, device=device)
        labels = torch.cat([torch.ones(len(bu)), torch.zeros(len(bu))]).to(device)

        opt.zero_grad()
        logits = model(u_t, v_t, r_t)
        loss = loss_fn(logits, labels)
        loss.backward()
        opt.step()
        t_epoch = time.time() - te0

        if epoch % check_every == 0 or epoch == 1:
            model.eval()
            with torch.no_grad():
                vscores = torch.sigmoid(model(val_u_t, val_v_t, val_r_t)).cpu().numpy()
                val_auprc = average_precision_score(val_labels, vscores)
            history.append({"epoch": epoch, "loss": float(loss.item()), "val_auprc": float(val_auprc),
                             "t_epoch": t_epoch})
            _log(f"epoch {epoch}/{max_epochs} loss={loss.item():.4f} val_auprc={val_auprc:.4f} "
                 f"t_epoch={t_epoch:.2f}s t_total={time.time()-t0:.1f}s")
            if val_auprc > best_val_auprc:
                best_val_auprc = val_auprc
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
                epochs_no_improve = 0
            else:
                epochs_no_improve += check_every
                if epochs_no_improve >= patience:
                    _log(f"Early stopping en epoch {epoch}. Mejor val_auprc={best_val_auprc:.4f}")
                    break

    model.load_state_dict(best_state)
    model.eval()
    z_final = model.node_emb.weight.detach().cpu().numpy()
    embeddings = {node_list[i]: z_final[i] for i in range(n_nodes)}
    timings = {
        "t_total": time.time() - t0,
        "n_epochs_run": history[-1]["epoch"] if history else 0,
        "best_val_auprc_interno": best_val_auprc,
        "metapaths_usados": metapaths_presentes,
        "arquitectura": {"dim": dim, "walk_length": walk_length, "lr": lr,
                          "weight_decay": weight_decay, "batch_size": batch_size},
        "history": history,
    }
    return embeddings, timings
