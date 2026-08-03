#!/usr/bin/env python3
"""
Seccion 4 (TFM graph_v4): GAT Autoencoder (GAT-AE), implementado a
mano en PyTorch puro (sin torch_geometric, misma decision que GAE).

IMPORTANTE (nomenclatura, exigido por el protocolo): este modelo usa
atencion (GAT) pero NO incorpora explicitamente tipos de nodo ni de
relacion en el mecanismo de atencion -- todas las aristas se tratan de
forma homogenea. Por lo tanto NO se lo denomina "heterogeneo" en
ningun documento de la memoria; es un autoencoder con atencion sobre un
grafo tratado como homogeneo, no una extension tipo HAN/HGT.

CARACTERISTICAS DE NODO:
    Condicion 'structural' (comportamiento original, Secciones 4-9):
    identico a GAE -- embedding entrenable por nodo, sin atributos
    biologicos.
    Condiciones 'uniprot_fixed' / 'uniprot_finetuned' (Seccion 10):
    identico a GAE -- ver NodeFeatureModule en gae_common.py (modulo
    compartido entre ambos encoders). Los nodos GEN con cobertura
    UniProt (80.18%, Seccion 8) inician desde su vector ProtT5 de
    1024-d, proyectado por una capa aprendible; el resto de nodos
    siguen con la tabla nn.Embedding entrenable de respaldo.

POR QUE EL SUBGRAFO COMUN (no el grafo de entrenamiento completo):
a diferencia de la multiplicacion dispersa de GAE (torch.sparse.mm,
que scipy/torch resuelven sin materializar un tensor por arista), el
mecanismo de atencion de GAT necesita materializar, para cada arista y
cada cabeza de atencion, un vector de mensaje de dimension out_dim.
Con el grafo de entrenamiento completo (~220M aristas simetrizadas) y
4 cabezas x 32 dim en la primera capa, eso son ~220M x 4 x 32 x 4
bytes =~ 112 GB solo para ese tensor intermedio -- inviable en la GPU
disponible (Tesla P100, 12GB). Con el subgrafo comun (~6.75M aristas
simetrizadas + auto-bucles) el mismo calculo da ~3.5GB, perfectamente
viable. Por eso GAT-AE entrena directamente sobre el subgrafo comun,
a diferencia de GAE (que si corrio sobre el grafo completo).

ARQUITECTURA:
  - Capa 1: GAT multi-cabeza (4 cabezas, 32 dim/cabeza, concatenadas
    -> 128 dim), activacion ELU, dropout 0.2.
  - Capa 2: GAT de 1 cabeza, 128 dim de salida (dimension final del
    embedding, sin concatenar -- una sola cabeza no la necesita).
  - Softmax de atencion implementado a mano via scatter (reduce='amax'
    + exp + scatter_add), sin depender de torch_geometric/torch_scatter.

ENTRENAMIENTO Y VALIDACION: identico protocolo que GAE (mismos splits,
misma dimension final, mismos negativos 1:1 muestreados por epoca,
misma validacion interna de reconstruccion con early stopping por
AUPRC), tal como exige el protocolo para GAT-AE.

LIMITE COMPUTACIONAL ADICIONAL (documentado, job 308498): incluso sobre
el subgrafo comun, la estimacion analitica de ~3.5GB fue optimista --
el tensor de mensajes se materializa en AMBAS capas (no solo la
primera), y ademas el grafo de autograd retiene el tensor de la capa 1
para poder calcular su gradiente durante el backward de la capa 2. Con
E=6,652,492 aristas (incl. auto-bucles), capa1 (4 cabezas x 32 dim) y
capa2 (1 cabeza x 128 dim) cada una materializa ~3.4GB en fp32, y la
suma mas el resto de tensores de activacion excede los 12GB del P100
(fallo real: 10.85GB en uso + 3.17GB solicitados > 11.90GB totales).
Fix aplicado: entrenamiento en precision mixta (torch.amp autocast a
float16 + GradScaler) en lugar de reducir la arquitectura -- esto
preserva exactamente la dimension 128 y el numero de cabezas exigidos
por el protocolo ("mismo protocolo que GAE"), a costa de una perdida de
precision numerica estandar en el forward/backward (los pesos maestros
y el optimizador siguen en fp32). Alternativa descartada: reducir
heads1/hidden_dim, porque eso cambiaria la arquitectura documentada sin
necesidad, cuando la causa real es memoria de activacion, no de
parametros.

SEGUNDO PROBLEMA DESCUBIERTO Y CORREGIDO (documentado, job 308501): la
version con precision mixta corrio sin error de memoria pero NO
entreno -- loss identica a ln(2)=0.6931 en todas las epocas registradas
(1, 5, 10) y early stopping a la epoca 10 con val_auprc~0.50 (nivel
azar), lo que en evaluacion de test colapso a auprc=0.0909/auroc=0.5000
(score constante) en 4 de 7 tareas. Diagnostico: GradScaler salta el
paso de optimizacion completo cuando detecta inf/nan en el gradiente
escalado; en un grafo con nodos hub (grado alto incluso en el subgrafo
comun), el gradiente de la capa de atencion se acumula via
index_add/scatter sobre potencialmente miles de aristas por nodo, y esa
suma, multiplicada por el factor de escala de fp16, puede superar el
rango representable -- si esto ocurre en TODAS las epocas, los pesos
nunca se actualizan desde su inicializacion aleatoria (consistente con
la perdida exactamente constante observada). Fix: gradient clipping
(clip_grad_norm_ a norma maxima 1.0) aplicado DESPUES de
scaler.unscale_() y ANTES de scaler.step(), que acota la norma del
gradiente antes de que el scaler decida si el paso es valido. Se
instrumento ademas el log por epoca con grad_norm, amp_scale y
step_skipped_amp para verificar empiricamente si el problema persiste.

TERCER PROBLEMA, DIAGNOSTICO DEFINITIVO (documentado, job 308504): con
grad_norm en notacion cientifica y patience=40 se confirmo que NO es un
bug de grafo desconectado -- el gradiente es real pero decae de forma
exponencial y monotonica (7.0e-07 en epoca 1 -> 9.7e-08 (5) -> 5.3e-09
(10) -> 1.2e-10 (15) -> ~5e-11 (20-30) -> 0.0 exacto (35-40)),
consistente con un colapso de representacion: el embedding entrenable
compartido, promediado dos veces por softmax de atencion (una operacion
inherentemente contractiva, sobre todo en nodos hub de grado muy alto),
se encoge rapidamente hacia el origen bajo la presion de weight_decay
(5e-4, igual que GAE) mientras la señal de gradiente que lo empuja hacia
afuera (proporcional a las normas de z, via el decodificador producto
punto) se encoge en la misma proporcion -- un punto fijo degenerado
autoreforzado en z=0, donde logit=0 para todo par (de ahi loss=ln(2)
exacto desde la epoca 1). GAE no sufre esto porque su propagacion GCN
(normalizacion simetrica fija D^-1/2(A+I)D^-1/2) no es una contraccion
tan agresiva como el promedio ponderado por atencion en dos capas.
Fix: weight_decay=0.0 especificamente para GAT-AE (no para GAE, que
funciona bien con 5e-4) -- elimina la presion de decaimiento que
domina a la señal de gradiente ya de por si debil al inicio, sin tocar
dimension/arquitectura/heads exigidos por el protocolo. Se agrega
tambien el log de la norma media de los embeddings (z_norm) por epoca
para verificar empiricamente que el colapso no vuelve a ocurrir.

Seccion 10 -- nota sobre uniprot_finetuned + weight_decay: el mismo
riesgo de colapso de representacion documentado arriba (parametro
compartido erosionado por weight_decay hasta z=0) motiva que, en modo
'finetuned', uniprot_raw se entrene en un grupo de parametros SIN
weight_decay (ver _build_optimizer en gae_common.py), incluso aunque
weight_decay ya sea 0.0 para el resto de GAT-AE por el fix anterior --
asi el comportamiento es correcto y explicito independientemente del
valor de weight_decay que reciba la funcion.
"""
import time

import numpy as np
import torch
import torch.nn as nn
from torch.amp import autocast, GradScaler
from sklearn.metrics import average_precision_score

from gae_common import build_graph, sample_negatives_numba, load_uniprot_features, NodeFeatureModule, _build_optimizer


def log(msg):
    print(f"[gat] {msg}", flush=True)


def build_edge_index(A, device):
    """edge_index (2,E) con auto-bucles, para message passing GAT."""
    coo = A.tocoo()
    n = A.shape[0]
    self_loops = np.arange(n)
    src = np.concatenate([coo.row, self_loops])
    dst = np.concatenate([coo.col, self_loops])
    edge_index = torch.tensor(np.stack([src, dst]), dtype=torch.long, device=device)
    return edge_index


def scatter_softmax(e, index, num_nodes):
    e_max = torch.full((num_nodes,), float("-inf"), device=e.device, dtype=e.dtype)
    e_max = e_max.scatter_reduce(0, index, e, reduce="amax", include_self=True)
    e_max_per_edge = e_max[index]
    e_exp = torch.exp(e - e_max_per_edge)
    denom = torch.zeros(num_nodes, device=e.device, dtype=e.dtype)
    denom = denom.scatter_add(0, index, e_exp)
    denom_per_edge = denom[index] + 1e-16
    return e_exp / denom_per_edge


class GATLayer(nn.Module):
    def __init__(self, in_dim, out_dim, heads=4, concat=True, dropout=0.2, negative_slope=0.2):
        super().__init__()
        self.heads = heads
        self.out_dim = out_dim
        self.concat = concat
        self.W = nn.Linear(in_dim, heads * out_dim, bias=False)
        self.a_src = nn.Parameter(torch.empty(heads, out_dim))
        self.a_dst = nn.Parameter(torch.empty(heads, out_dim))
        nn.init.xavier_uniform_(self.a_src)
        nn.init.xavier_uniform_(self.a_dst)
        self.leakyrelu = nn.LeakyReLU(negative_slope)
        self.dropout = nn.Dropout(dropout)

    def forward(self, H, edge_index, num_nodes):
        src, dst = edge_index[0], edge_index[1]
        Wh = self.W(H).view(-1, self.heads, self.out_dim)  # (N, heads, out_dim)

        alpha_src = (Wh * self.a_src).sum(dim=-1)  # (N, heads)
        alpha_dst = (Wh * self.a_dst).sum(dim=-1)  # (N, heads)
        e = self.leakyrelu(alpha_src[src] + alpha_dst[dst])  # (E, heads)

        alpha = torch.empty_like(e)
        for h in range(self.heads):
            alpha[:, h] = scatter_softmax(e[:, h].contiguous(), dst, num_nodes)
        alpha = self.dropout(alpha)
        # softmax se calcula en la precision de e (fp32, estable) pero se convierte
        # explicitamente al dtype de Wh (fp16 bajo autocast) ANTES de multiplicar por
        # Wh[src]: si no, la promocion automatica fp16*fp32 -> fp32 de PyTorch anularia
        # el ahorro de memoria buscado en el tensor grande (E, heads, out_dim).
        alpha = alpha.to(Wh.dtype)

        messages = Wh[src] * alpha.unsqueeze(-1)  # (E, heads, out_dim)
        # dtype = messages.dtype (no H.dtype): bajo autocast, Wh/messages quedan en
        # fp16 aunque H (p.ej. el embedding de entrada, no afectado por autocast) siga
        # en fp32 -- index_add exige que out y messages compartan dtype.
        out = torch.zeros(num_nodes, self.heads, self.out_dim, device=H.device, dtype=messages.dtype)
        out = out.index_add(0, dst, messages)

        if self.concat:
            out = out.reshape(num_nodes, self.heads * self.out_dim)
        else:
            out = out.mean(dim=1)
        return out


class GATEncoder(nn.Module):
    """
    Features de entrada: ver NodeFeatureModule (gae_common.py, modulo
    compartido con GAE).
    """
    def __init__(self, n_nodes, in_dim=128, hidden_dim=32, heads1=4, out_dim=128, dropout=0.2,
                 uniprot_vectors=None, uniprot_node_indices=None, uniprot_mode=None):
        super().__init__()
        self.feature_module = NodeFeatureModule(n_nodes, in_dim, uniprot_vectors,
                                                 uniprot_node_indices, uniprot_mode)
        self.gat1 = GATLayer(in_dim, hidden_dim, heads=heads1, concat=True, dropout=dropout)
        self.gat2 = GATLayer(hidden_dim * heads1, out_dim, heads=1, concat=False, dropout=dropout)
        self.elu = nn.ELU()

    def forward(self, edge_index, num_nodes):
        # self.feature_module() se ejecuta DENTRO del bloque autocast del
        # llamador (ver train_gat): la proyeccion UniProt (nn.Linear) queda
        # en fp16 igual que el resto del forward, y el cast explicito a
        # H.dtype dentro de NodeFeatureModule.forward evita el mismo tipo
        # de error de dtype ya resuelto en GATLayer (alpha.to(Wh.dtype)).
        H = self.feature_module()
        H = self.gat1(H, edge_index, num_nodes)
        H = self.elu(H)
        H = self.gat2(H, edge_index, num_nodes)
        return H


def decode(z, u_idx, v_idx):
    return (z[u_idx] * z[v_idx]).sum(dim=1)


def train_gat(train_graph_path, dim=128, heads1=4, hidden_dim=32, dropout=0.2,
              lr=0.005, weight_decay=5e-4, max_epochs=200, patience=10, check_every=5,
              batch_edges=200_000, n_val_internal=20_000, seed=42, device=None, log_prefix="gat",
              use_amp=True, grad_clip_norm=1.0, uniprot_tsv=None, uniprot_mode=None):
    def _log(msg):
        print(f"[{log_prefix}] {msg}", flush=True)

    torch.manual_seed(seed)
    np.random.seed(seed)
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Device: {device}")
    amp_enabled = use_amp and device.type == "cuda"
    _log(f"Precision mixta (autocast fp16 + GradScaler): {'activada' if amp_enabled else 'desactivada'}")
    scaler = GradScaler(device.type if device.type == "cuda" else "cpu", enabled=amp_enabled)

    t0 = time.time()
    A, pos_pairs, nodes, idx = build_graph(train_graph_path)
    n_nodes = len(nodes)
    _log(f"{n_nodes:,} nodos, {A.nnz:,} no-ceros (simetrizado), "
         f"{len(pos_pairs):,} aristas originales (una direccion). t={time.time() - t0:.1f}s")

    uniprot_vectors, uniprot_node_indices = None, None
    if uniprot_tsv is not None:
        _log(f"Cargando features externas UniProt (modo={uniprot_mode}) desde {uniprot_tsv}...")
        uniprot_node_indices, uniprot_vectors = load_uniprot_features(uniprot_tsv, idx)

    edge_index = build_edge_index(A, device)
    _log(f"edge_index construido ({edge_index.shape[1]:,} aristas incl. auto-bucles). t={time.time() - t0:.1f}s")

    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(pos_pairs))
    val_idx_internal = perm[:n_val_internal]
    train_pool_idx = perm[n_val_internal:]
    val_pos = pos_pairs[val_idx_internal]
    _log(f"Validacion interna de reconstruccion (GAT-AE): {len(val_pos):,} aristas positivas")

    val_neg_u, val_neg_v = sample_negatives_numba(A.indptr.astype(np.int64), A.indices.astype(np.int64),
                                                    n_nodes, len(val_pos), seed + 999)

    model = GATEncoder(n_nodes, in_dim=dim, hidden_dim=hidden_dim, heads1=heads1,
                        out_dim=dim, dropout=dropout,
                        uniprot_vectors=uniprot_vectors, uniprot_node_indices=uniprot_node_indices,
                        uniprot_mode=uniprot_mode).to(device)
    initial_projection = model.feature_module.initial_projection_snapshot()
    opt = _build_optimizer(model, lr, weight_decay, uniprot_mode, _log)
    loss_fn = nn.BCEWithLogitsLoss()

    best_val_auprc = -1.0
    best_state = None
    epochs_no_improve = 0
    history = []

    for epoch in range(1, max_epochs + 1):
        te0 = time.time()
        model.train()
        opt.zero_grad()

        batch_pos_idx = rng.choice(train_pool_idx, size=min(batch_edges, len(train_pool_idx)), replace=False)
        batch_pos = pos_pairs[batch_pos_idx]
        neg_u, neg_v = sample_negatives_numba(A.indptr.astype(np.int64), A.indices.astype(np.int64),
                                               n_nodes, len(batch_pos), seed + epoch)

        u_idx = torch.tensor(np.concatenate([batch_pos[:, 0], neg_u]), dtype=torch.long, device=device)
        v_idx = torch.tensor(np.concatenate([batch_pos[:, 1], neg_v]), dtype=torch.long, device=device)
        labels = torch.cat([torch.ones(len(batch_pos)), torch.zeros(len(neg_u))]).to(device)

        with autocast(device.type, dtype=torch.float16, enabled=amp_enabled):
            z = model(edge_index, n_nodes)
            logits = decode(z, u_idx, v_idx)
            loss = loss_fn(logits.float(), labels)
        scale_before = scaler.get_scale()
        scaler.scale(loss).backward()
        # Sin esto, un solo nodo hub (grado muy alto incluso en el subgrafo
        # comun) puede acumular, via index_add/scatter, un gradiente cuya
        # version escalada (para fp16) supera el rango representable y se
        # vuelve inf/nan -- GradScaler entonces salta el paso por completo.
        # Si esto ocurre en TODAS las epocas el modelo nunca se actualiza
        # (perdida identica a ln(2) para siempre, visto empiricamente en el
        # job 308501 sin este clipping). unscale_+clip_grad_norm_ acota la
        # norma del gradiente ANTES de que GradScaler decida si el paso es
        # valido, evitando ese colapso sin perder la aceleracion de fp16.
        scaler.unscale_(opt)
        grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
        scaler.step(opt)
        scaler.update()
        step_skipped = scaler.get_scale() < scale_before
        t_epoch = time.time() - te0

        if epoch % check_every == 0 or epoch == 1:
            model.eval()
            with torch.no_grad(), autocast(device.type, dtype=torch.float16, enabled=amp_enabled):
                z_eval = model(edge_index, n_nodes)
                vu = torch.tensor(np.concatenate([val_pos[:, 0], val_neg_u]), dtype=torch.long, device=device)
                vv = torch.tensor(np.concatenate([val_pos[:, 1], val_neg_v]), dtype=torch.long, device=device)
                vlabels = np.concatenate([np.ones(len(val_pos)), np.zeros(len(val_neg_u))])
                vscores = torch.sigmoid(decode(z_eval, vu, vv).float()).cpu().numpy()
                val_auprc = average_precision_score(vlabels, vscores)
                z_norm = z_eval.detach().float().norm(dim=1).mean().item()

            history.append({"epoch": epoch, "loss": float(loss.item()), "val_auprc": float(val_auprc),
                             "t_epoch": t_epoch, "grad_norm": float(grad_norm), "z_norm": float(z_norm),
                             "amp_scale": float(scaler.get_scale()), "step_skipped_amp": bool(step_skipped)})
            _log(f"epoch {epoch}/{max_epochs} loss={loss.item():.6f} val_auprc={val_auprc:.4f} "
                 f"grad_norm={float(grad_norm):.3e} z_norm={z_norm:.4e} amp_scale={scaler.get_scale():.0f} "
                 f"step_skipped={step_skipped} t_epoch={t_epoch:.2f}s t_total={time.time() - t0:.1f}s")

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
    with torch.no_grad():
        z_final = model(edge_index, n_nodes).cpu().numpy()

    embeddings = {nodes[i]: z_final[i] for i in range(n_nodes)}

    initial_uniprot_embeddings = None
    if initial_projection is not None:
        initial_uniprot_embeddings = {nodes[i]: initial_projection[j]
                                       for j, i in enumerate(uniprot_node_indices)}

    timings = {
        "t_total": time.time() - t0,
        "n_epochs_run": history[-1]["epoch"] if history else 0,
        "best_val_auprc_interno": best_val_auprc,
        "arquitectura": {
            "capa1_heads": heads1, "capa1_dim_por_cabeza": hidden_dim,
            "capa2_heads": 1, "capa2_dim_salida": dim,
            "dropout": dropout, "lr": lr, "weight_decay": weight_decay,
            "precision_mixta_fp16": amp_enabled, "grad_clip_norm": grad_clip_norm,
        },
        "history": history,
        "uniprot_mode": uniprot_mode,
        "n_uniprot_covered": int(len(uniprot_node_indices)) if uniprot_node_indices is not None else 0,
        "initial_uniprot_embeddings": initial_uniprot_embeddings,
    }
    return embeddings, timings
