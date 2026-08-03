#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): GAE (Graph Autoencoder) implementado a mano
en PyTorch puro, sin torch_geometric (decision confirmada con el
usuario: evitar la instalacion de una dependencia nueva cuya
disponibilidad de internet en los nodos de computo no esta garantizada).

CARACTERISTICAS DE NODO USADAS (documentacion explicita, exigida por
el protocolo -- "no inventes atributos biologicos"):
    Condicion 'structural' (comportamiento original, Secciones 4-9):
    NO se usan atributos biologicos reales. Se usa un EMBEDDING
    ENTRENABLE por nodo (nn.Embedding inicializado aleatoriamente) como
    "features" de entrada al encoder -- el enfoque estandar de GAE
    (Kipf & Welling 2016) cuando no hay atributos de nodo reales
    disponibles.
    Condiciones 'uniprot_fixed' / 'uniprot_finetuned' (Seccion 10):
    los nodos GEN con cobertura UniProt (21.111/26.330 = 80.18%, ver
    Seccion 8) inician desde su vector ProtT5 de 1024-d real, proyectado
    mediante una capa nn.Linear(1024, dim) aprendible (NodeFeatureModule,
    definido en este mismo archivo y reutilizado por gat_common.py). El
    resto de nodos -- no-GEN y GEN SIN cobertura UniProt -- siguen
    usando la tabla nn.Embedding entrenable de respaldo: nunca se
    inventa un vector externo para ellos (regla explicita, Seccion 8,
    Seccion 9 punto 5). 'fixed' congela la tabla UniProt cruda
    (requires_grad=False); 'finetuned' permite que tambien se ajuste
    durante el entrenamiento.

ARQUITECTURA: 2 capas GCN dispersas (D^-1/2 (A+I) D^-1/2 @ H @ W),
implementadas con torch.sparse.mm (nunca se materializa una adyacencia
densa). Decoder de producto interno (score(u,v) = z_u . z_v), estandar
de GAE.

ENTRENAMIENTO: la propagacion (forward pass del encoder) usa TODO el
grafo de aristas de entrenamiento (message passing completo); la
perdida de reconstruccion (BCE con logits) se calcula sobre un
subconjunto muestreado de aristas positivas + negativos muestreados en
cada epoca (no las ~110M aristas completas en cada paso, por costo).
Se aplica perdida ponderada (pos_weight en BCEWithLogitsLoss) si hay
desbalance entre el conteo de positivos y negativos de la epoca (aqui
1:1, por lo que pos_weight=1.0 en la config por defecto, pero queda
parametrizado por si se cambia el ratio).

VALIDACION INTERNA / EARLY STOPPING: se reserva una muestra pequena y
fija de aristas de entrenamiento (nunca las aristas de val/test de las
7 tareas primarias, que ya estan excluidas de graph_edges_train_only.tsv
desde la Seccion 3) como validacion INTERNA de la tarea de
reconstruccion del propio GAE, para AUPRC de validacion / early
stopping, tal como exige el protocolo. Esto es distinto de la
evaluacion final en las 7 tareas primarias (que se hace aparte, con el
mismo common_eval.py que usan los demas metodos).
"""
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from numba import njit, prange
from scipy.sparse import csr_matrix


def log(msg):
    print(f"[gae] {msg}", flush=True)


def build_graph(train_graph_path):
    df = pd.read_csv(train_graph_path, sep="\t", usecols=["node1", "node2"], dtype=str)
    nodes = sorted(set(df["node1"]) | set(df["node2"]))
    idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)

    rows = df["node1"].map(idx).to_numpy()
    cols = df["node2"].map(idx).to_numpy()
    all_rows = np.concatenate([rows, cols]).astype(np.int64)
    all_cols = np.concatenate([cols, rows]).astype(np.int64)
    data = np.ones(len(all_rows), dtype=np.int8)
    A = csr_matrix((data, (all_rows, all_cols)), shape=(n_nodes, n_nodes))
    A.data[:] = 1
    A.sum_duplicates()
    A.data[:] = 1
    A.sort_indices()

    pos_pairs = np.stack([rows, cols], axis=1)  # aristas originales (una direccion), para supervision
    return A, pos_pairs, nodes, idx


def load_uniprot_features(uniprot_tsv_path, idx):
    """
    Seccion 10: carga los embeddings crudos de UniProt (Seccion 8,
    formato node_id|node_type|embedding_source|dimension|e_0001...e_n)
    y los alinea con el indice de nodos (idx) de ESTE grafo/subgrafo.

    Devuelve:
      uniprot_node_indices: array int64, posiciones (segun idx) de los
        nodos cubiertos, en el mismo orden que uniprot_vectors.
      uniprot_vectors: array float32 (k, dim_externa).

    Los nodos GEN sin cobertura (no aparecen en este archivo) y todos
    los nodos no-GEN quedan fuera deliberadamente -- nunca se inventa
    un vector externo para ellos (regla explicita de la Seccion 8);
    siguen usando la tabla nn.Embedding entrenable de respaldo en
    NodeFeatureModule. La cobertura especifica de GEN (80.18%) ya esta
    auditada y documentada en la Seccion 8; aqui solo se reporta la
    cobertura relativa al grafo/subgrafo actual, para no duplicar ese
    calculo con una heuristica fragil de tipo de nodo.
    """
    df = pd.read_csv(uniprot_tsv_path, sep="\t")
    e_cols = sorted((c for c in df.columns if c.startswith("e_")), key=lambda c: int(c.split("_")[1]))
    node_indices, vectors, n_missing = [], [], 0
    vecs_arr = df[e_cols].to_numpy(dtype=np.float32)
    for node_id, vec in zip(df["node_id"].tolist(), vecs_arr):
        if node_id in idx:
            node_indices.append(idx[node_id])
            vectors.append(vec)
        else:
            n_missing += 1
    if n_missing:
        log(f"AVISO: {n_missing:,} nodos con embedding UniProt no estan en este grafo (se ignoran).")
    uniprot_node_indices = np.asarray(node_indices, dtype=np.int64)
    uniprot_vectors = (np.stack(vectors).astype(np.float32) if vectors
                       else np.zeros((0, len(e_cols)), dtype=np.float32))
    log(f"UniProt cargado y alineado: {len(uniprot_node_indices):,} nodos cubiertos de "
        f"{len(idx):,} nodos totales del grafo ({100 * len(uniprot_node_indices) / len(idx):.2f}%). "
        f"Dim externa={len(e_cols)}. (Cobertura especifica de GEN: 21.111/26.330=80.18%, Seccion 8.)")
    return uniprot_node_indices, uniprot_vectors


def normalized_adj_torch(A, device):
    """D^-1/2 (A+I) D^-1/2 como torch.sparse_coo_tensor."""
    n = A.shape[0]
    A_self = A + csr_matrix((np.ones(n, dtype=np.float32), (np.arange(n), np.arange(n))), shape=(n, n))
    A_self = A_self.tocoo()
    deg = np.asarray(A_self.sum(axis=1)).ravel()
    deg_inv_sqrt = np.power(deg, -0.5, where=deg > 0)
    deg_inv_sqrt[deg == 0] = 0.0
    vals = deg_inv_sqrt[A_self.row] * A_self.data * deg_inv_sqrt[A_self.col]

    indices = torch.tensor(np.stack([A_self.row, A_self.col]), dtype=torch.long)
    values = torch.tensor(vals, dtype=torch.float32)
    A_norm = torch.sparse_coo_tensor(indices, values, (n, n)).coalesce().to(device)
    return A_norm


@njit(cache=True)
def _has_edge(indices, start, end, target):
    lo, hi = start, end
    while lo < hi:
        mid = (lo + hi) // 2
        v = indices[mid]
        if v == target:
            return True
        elif v < target:
            lo = mid + 1
        else:
            hi = mid
    return False


@njit(cache=True, parallel=True)
def sample_negatives_numba(indptr, indices, n_nodes, n_samples, seed):
    out_u = np.empty(n_samples, dtype=np.int64)
    out_v = np.empty(n_samples, dtype=np.int64)
    for i in prange(n_samples):
        np.random.seed(seed + i)
        while True:
            u = np.random.randint(0, n_nodes)
            v = np.random.randint(0, n_nodes)
            if u == v:
                continue
            start, end = indptr[u], indptr[u + 1]
            if _has_edge(indices, start, end, v):
                continue
            out_u[i] = u
            out_v[i] = v
            break
    return out_u, out_v


class SparseGCNLayer(nn.Module):
    def __init__(self, in_dim, out_dim, activation=True):
        super().__init__()
        self.lin = nn.Linear(in_dim, out_dim)
        self.activation = activation

    def forward(self, A_norm, H):
        H = self.lin(H)
        H = torch.sparse.mm(A_norm, H)
        if self.activation:
            H = torch.relu(H)
        return H


class NodeFeatureModule(nn.Module):
    """
    Seccion 10: genera la matriz de features de entrada (n_nodes, in_dim)
    para el encoder (GAE o GAT-AE). Compartido por gae_common.py y
    gat_common.py para no duplicar la logica entre ambos.

    Condicion 'structural' (uniprot_vectors=None): comportamiento
    identico al original -- una tabla nn.Embedding(n_nodes, in_dim)
    entrenable, init Xavier, sin atributos biologicos.

    Condiciones 'uniprot_fixed' / 'uniprot_finetuned': los nodos GEN con
    cobertura UniProt (ver load_uniprot_features) inician desde su
    vector UniProt crudo, pasado por una capa nn.Linear(dim_ext, in_dim)
    aprendible (gene_projection). El resto de nodos -- no-GEN y GEN SIN
    cobertura UniProt -- siguen usando la tabla nn.Embedding entrenable
    de respaldo: nunca se inventa un vector externo para ellos.

    'fixed' vs 'finetuned' controla si la tabla UniProt cruda
    (uniprot_raw) tambien recibe gradiente:
      - fixed:     buffer (requires_grad=False); solo gene_projection y
                   el resto del encoder aprenden.
      - finetuned: nn.Parameter (requires_grad=True); uniprot_raw
                   tambien se ajusta durante el entrenamiento.

    La combinacion con la tabla de respaldo usa index_copy (out-of-place),
    diferenciable respecto a ambas fuentes sin operaciones in-place que
    rompan el grafo de autograd. El cast explicito `.to(H.dtype)` antes
    de index_copy replica la misma precaucion de dtype ya necesaria en
    GATLayer bajo precision mixta (ver gat_common.py) -- aqui H siempre
    es fp32 (parametro de nn.Embedding, no afectado por autocast), pero
    `projected` si puede quedar en fp16 bajo autocast en GAT-AE.
    """

    def __init__(self, n_nodes, in_dim, uniprot_vectors=None, uniprot_node_indices=None,
                 uniprot_mode=None):
        super().__init__()
        self.node_features = nn.Embedding(n_nodes, in_dim)
        nn.init.xavier_uniform_(self.node_features.weight)
        self.has_uniprot = uniprot_vectors is not None and len(uniprot_vectors) > 0
        self.uniprot_mode = uniprot_mode
        if self.has_uniprot:
            assert uniprot_mode in ("fixed", "finetuned"), \
                "uniprot_mode debe ser 'fixed' o 'finetuned' si se proveen uniprot_vectors"
            ext_dim = uniprot_vectors.shape[1]
            self.gene_projection = nn.Linear(ext_dim, in_dim)
            nn.init.xavier_uniform_(self.gene_projection.weight)
            nn.init.zeros_(self.gene_projection.bias)
            uniprot_t = torch.as_tensor(uniprot_vectors, dtype=torch.float32)
            if uniprot_mode == "fixed":
                self.register_buffer("uniprot_raw", uniprot_t)
            else:  # finetuned
                self.uniprot_raw = nn.Parameter(uniprot_t.clone())
            self.register_buffer("uniprot_node_indices",
                                  torch.as_tensor(uniprot_node_indices, dtype=torch.long))

    def initial_projection_snapshot(self):
        """Vector proyectado (externo -> in_dim) ANTES de cualquier paso de
        entrenamiento -- 'embeddings iniciales' exigidos junto a los finales
        por el protocolo de la Seccion 10. None si no hay condicion UniProt."""
        if not self.has_uniprot:
            return None
        with torch.no_grad():
            return self.gene_projection(self.uniprot_raw).detach().clone().cpu().numpy()

    def forward(self):
        H = self.node_features.weight
        if self.has_uniprot:
            projected = self.gene_projection(self.uniprot_raw)
            projected = projected.to(H.dtype)
            H = H.index_copy(0, self.uniprot_node_indices, projected)
        return H


class GAEEncoder(nn.Module):
    """
    Features de entrada: ver NodeFeatureModule. 2 capas GCN dispersas.
    """
    def __init__(self, n_nodes, in_dim=128, hidden_dim=128, out_dim=128,
                 uniprot_vectors=None, uniprot_node_indices=None, uniprot_mode=None):
        super().__init__()
        self.feature_module = NodeFeatureModule(n_nodes, in_dim, uniprot_vectors,
                                                 uniprot_node_indices, uniprot_mode)
        self.gcn1 = SparseGCNLayer(in_dim, hidden_dim, activation=True)
        self.gcn2 = SparseGCNLayer(hidden_dim, out_dim, activation=False)

    def forward(self, A_norm):
        H = self.feature_module()
        H = self.gcn1(A_norm, H)
        H = self.gcn2(A_norm, H)
        return H


def decode(z, u_idx, v_idx):
    return (z[u_idx] * z[v_idx]).sum(dim=1)


def _build_optimizer(model, lr, weight_decay, uniprot_mode, _log):
    """
    Seccion 10: en modo 'finetuned', uniprot_raw se entrena SIN
    weight_decay (grupo de parametros separado) -- mismo razonamiento
    que ya se aplico a todo GAT-AE (ver gat_common.py): no queremos que
    el decaimiento erosione una senal biologica pre-entrenada real, solo
    queremos regularizar las partes inicializadas aleatoriamente.
    """
    if uniprot_mode == "finetuned":
        uniprot_params = [model.feature_module.uniprot_raw]
        uniprot_param_ids = {id(p) for p in uniprot_params}
        other_params = [p for p in model.parameters() if id(p) not in uniprot_param_ids]
        _log("Optimizador con grupos separados: weight_decay=0.0 para uniprot_raw "
             f"(preserva la senal biologica pre-entrenada), weight_decay={weight_decay} para el resto.")
        return torch.optim.Adam([
            {"params": other_params, "weight_decay": weight_decay},
            {"params": uniprot_params, "weight_decay": 0.0},
        ], lr=lr)
    return torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)


def train_gae(train_graph_path, dim=128, hidden_dim=128, lr=0.01, weight_decay=5e-4,
              max_epochs=200, patience=10, check_every=5, batch_edges=500_000,
              n_val_internal=50_000, seed=42, device=None, log_prefix="gae",
              uniprot_tsv=None, uniprot_mode=None):
    def _log(msg):
        print(f"[{log_prefix}] {msg}", flush=True)

    torch.manual_seed(seed)
    np.random.seed(seed)
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Device: {device}")

    t0 = time.time()
    A, pos_pairs, nodes, idx = build_graph(train_graph_path)
    n_nodes = len(nodes)
    _log(f"{n_nodes:,} nodos, {A.nnz:,} no-ceros (simetrizado), "
         f"{len(pos_pairs):,} aristas originales (una direccion). t={time.time() - t0:.1f}s")

    uniprot_vectors, uniprot_node_indices = None, None
    if uniprot_tsv is not None:
        _log(f"Cargando features externas UniProt (modo={uniprot_mode}) desde {uniprot_tsv}...")
        uniprot_node_indices, uniprot_vectors = load_uniprot_features(uniprot_tsv, idx)

    A_norm = normalized_adj_torch(A, device)
    _log(f"Adyacencia normalizada en torch (sparse) lista. t={time.time() - t0:.1f}s")

    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(pos_pairs))
    val_idx_internal = perm[:n_val_internal]
    train_pool_idx = perm[n_val_internal:]
    val_pos = pos_pairs[val_idx_internal]
    _log(f"Validacion interna de reconstruccion (GAE): {len(val_pos):,} aristas positivas "
         f"(excluidas SOLO de la perdida, no de la propagacion)")

    val_neg_u, val_neg_v = sample_negatives_numba(A.indptr.astype(np.int64), A.indices.astype(np.int64),
                                                    n_nodes, len(val_pos), seed + 999)

    model = GAEEncoder(n_nodes, in_dim=dim, hidden_dim=hidden_dim, out_dim=dim,
                        uniprot_vectors=uniprot_vectors, uniprot_node_indices=uniprot_node_indices,
                        uniprot_mode=uniprot_mode).to(device)
    initial_projection = model.feature_module.initial_projection_snapshot()
    opt = _build_optimizer(model, lr, weight_decay, uniprot_mode, _log)
    loss_fn = nn.BCEWithLogitsLoss()

    best_val_auprc = -1.0
    best_state = None
    epochs_no_improve = 0
    from sklearn.metrics import average_precision_score

    history = []
    for epoch in range(1, max_epochs + 1):
        te0 = time.time()
        model.train()
        opt.zero_grad()
        z = model(A_norm)

        batch_pos_idx = rng.choice(train_pool_idx, size=min(batch_edges, len(train_pool_idx)), replace=False)
        batch_pos = pos_pairs[batch_pos_idx]
        neg_u, neg_v = sample_negatives_numba(A.indptr.astype(np.int64), A.indices.astype(np.int64),
                                               n_nodes, len(batch_pos), seed + epoch)

        u_idx = torch.tensor(np.concatenate([batch_pos[:, 0], neg_u]), dtype=torch.long, device=device)
        v_idx = torch.tensor(np.concatenate([batch_pos[:, 1], neg_v]), dtype=torch.long, device=device)
        labels = torch.cat([torch.ones(len(batch_pos)), torch.zeros(len(neg_u))]).to(device)

        logits = decode(z, u_idx, v_idx)
        loss = loss_fn(logits, labels)
        loss.backward()
        opt.step()

        t_epoch = time.time() - te0

        if epoch % check_every == 0 or epoch == 1:
            model.eval()
            with torch.no_grad():
                z_eval = model(A_norm)
                vu = torch.tensor(np.concatenate([val_pos[:, 0], val_neg_u]), dtype=torch.long, device=device)
                vv = torch.tensor(np.concatenate([val_pos[:, 1], val_neg_v]), dtype=torch.long, device=device)
                vlabels = np.concatenate([np.ones(len(val_pos)), np.zeros(len(val_neg_u))])
                vscores = torch.sigmoid(decode(z_eval, vu, vv)).cpu().numpy()
                val_auprc = average_precision_score(vlabels, vscores)

            history.append({"epoch": epoch, "loss": float(loss.item()), "val_auprc": float(val_auprc),
                             "t_epoch": t_epoch})
            _log(f"epoch {epoch}/{max_epochs} loss={loss.item():.4f} val_auprc={val_auprc:.4f} "
                 f"t_epoch={t_epoch:.2f}s t_total={time.time() - t0:.1f}s")

            if val_auprc > best_val_auprc:
                best_val_auprc = val_auprc
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
                epochs_no_improve = 0
            else:
                epochs_no_improve += check_every
                if epochs_no_improve >= patience:
                    _log(f"Early stopping en epoch {epoch} (sin mejora en {epochs_no_improve} epocas). "
                         f"Mejor val_auprc={best_val_auprc:.4f}")
                    break

    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        z_final = model(A_norm).cpu().numpy()

    embeddings = {nodes[i]: z_final[i] for i in range(n_nodes)}

    initial_uniprot_embeddings = None
    if initial_projection is not None:
        initial_uniprot_embeddings = {nodes[i]: initial_projection[j]
                                       for j, i in enumerate(uniprot_node_indices)}

    timings = {"t_total": time.time() - t0,
               "n_epochs_run": history[-1]["epoch"] if history else 0,
               "best_val_auprc_interno": best_val_auprc, "history": history,
               "uniprot_mode": uniprot_mode,
               "n_uniprot_covered": int(len(uniprot_node_indices)) if uniprot_node_indices is not None else 0,
               "initial_uniprot_embeddings": initial_uniprot_embeddings}
    return embeddings, timings
