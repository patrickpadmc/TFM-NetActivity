"""
walker.py -- reescrito para TFM-NetActivity graph_v3, ya no basado en
e-yi/hin2vec_pytorch tal cual (aunque conserva su interfaz publica).

Motivo del reescrito: la version original (y el parche anterior con
MultiDiGraph de networkx) usa una estructura de diccionarios anidados
por arista, con overhead de varios cientos de bytes a >1KB por arista
en Python puro, y ademas guarda la adyacencia dos veces (sucesores y
predecesores). Con graph_v3 (varios millones de aristas incluso
despues de filtrar/submuestrear) esto satura la memoria del job antes
de llegar a entrenar.

Esta version representa el grafo como arrays de numpy en formato CSR
(offsets + vecinos), igual que build_adjacency() en spectral_embedding_v3.py:
~12 bytes por arista en vez de >1KB. model.py y train_hin2vec_v3.py no
necesitan cambios: la interfaz publica (window, path_size, id2node,
id2path, sample()) es la misma.

Los metapaths se siguen indexando por la secuencia de RELACIONES
(edge_class) recorridas en el walk, no por tipos de nodo -- ver
justificacion en versiones anteriores de este archivo / la
conversacion del proyecto.
"""
import numpy as np
import pandas as pd
from collections import defaultdict


class HIN:
    """
    Representa el grafo como CSR (offsets + vecinos) y genera vertex
    sequences (random walks) sobre esa estructura.
    """

    def __init__(self, window=None):
        self._window = window
        self._path2id = None
        self._id2path = None

        self.node_size = 0
        self._id2node = None
        self._id2edgeclass = None

        # CSR: para el nodo u, sus aristas salientes son
        # dst[indptr[u]:indptr[u+1]], edge_class_id[...], weight[...]
        self._indptr = None
        self._dst = None
        self._edge_class_id = None
        self._weight = None

    @property
    def id2node(self):
        return self._id2node

    @property
    def id2path(self):
        return self._id2path

    @property
    def window(self):
        return self._window

    @window.setter
    def window(self, val):
        if not self._window:
            self._window = val
        else:
            raise ValueError("window solo puede ser asignado una vez")

    @property
    def path_size(self):
        if self._path2id is None:
            raise ValueError("correr sample() primero para contar path size")
        return len(self._path2id)

    def build_from_edges(self, edges_df):
        """
        edges_df: DataFrame con columnas source_node, source_class,
        dest_node, dest_class, edge_class, weight. Se asume que YA
        incluye ambas direcciones (load_a_HIN_from_pandas se encarga
        de eso antes de llamar a este metodo).
        """
        print("factorizando nodos...", flush=True)
        combined = np.concatenate([
            edges_df["source_node"].to_numpy(),
            edges_df["dest_node"].to_numpy(),
        ])
        codes, uniques = pd.factorize(combined, sort=False)
        n = len(edges_df)
        src = codes[:n].astype(np.int32)
        dst = codes[n:].astype(np.int32)
        self.node_size = len(uniques)
        self._id2node = dict(enumerate(uniques))
        print(f"{self.node_size} nodos unicos", flush=True)

        print("factorizando edge_class (relaciones)...", flush=True)
        ec_codes, ec_uniques = pd.factorize(edges_df["edge_class"].to_numpy(), sort=False)
        self._id2edgeclass = dict(enumerate(ec_uniques))
        edge_class_id = ec_codes.astype(np.int32)
        print(f"{len(ec_uniques)} edge_class distintos", flush=True)

        weight = edges_df["weight"].to_numpy(dtype=np.float32)

        print("construyendo CSR...", flush=True)
        order = np.argsort(src, kind="stable")
        src_sorted = src[order]
        self._dst = dst[order]
        self._edge_class_id = edge_class_id[order]
        self._weight = weight[order]
        self._indptr = np.searchsorted(src_sorted, np.arange(self.node_size + 1))
        print(f"CSR listo: {len(self._dst)} aristas dirigidas", flush=True)

    def small_walk(self, start_node, length, rng):
        """
        Genera un random walk desde start_node.

        Returns:
            node_walk: lista de node ids visitados
            relation_walk: lista de edge_class (str) recorridos, tal que
                relation_walk[k] es la relacion usada para ir de
                node_walk[k] a node_walk[k+1].
        """
        node_walk = [start_node]
        relation_walk = []
        cur = start_node
        for _ in range(1, length):
            lo, hi = self._indptr[cur], self._indptr[cur + 1]
            if hi <= lo:
                break
            w = self._weight[lo:hi]
            probs = w / w.sum()
            idx = lo + rng.choice(hi - lo, p=probs)
            next_node = int(self._dst[idx])
            rel = self._id2edgeclass[int(self._edge_class_id[idx])]
            node_walk.append(next_node)
            relation_walk.append(rel)
            cur = next_node
        return node_walk, relation_walk

    def do_walks(self, length, rng):
        for start_node in range(self.node_size):
            yield self.small_walk(start_node, length, rng)

    def sample(self, length, n_repeat, max_samples=None, seed=None):
        """
        Genera muestras (start_id, end_id, path_id) recorriendo cada
        nodo como punto de partida n_repeat veces.

        Si max_samples se especifica, usa reservoir sampling (algoritmo R)
        para mantener una muestra uniforme de tamano max_samples SIN
        materializar nunca la lista completa de muestras crudas en
        memoria. El vocabulario de paths (self._path2id) se construye
        sobre el total de muestras vistas, no solo sobre las conservadas.

        :param length: longitud del walk
        :param n_repeat: numero de walks por nodo de partida
        :param max_samples: tamano maximo del reservoir (None = sin limite)
        :param seed: semilla, para reproducibilidad (walk + reservoir)
        :return: lista de (start_id, end_id, path_id)
        """
        if not self.window:
            raise ValueError("window no asignado")

        if self._path2id is None:
            self._path2id = defaultdict(lambda: len(self._path2id))

        rng = np.random.default_rng(seed)
        reservoir = []
        n_seen = 0

        def emit(item):
            nonlocal n_seen
            n_seen += 1
            if max_samples is None:
                reservoir.append(item)
            elif len(reservoir) < max_samples:
                reservoir.append(item)
            else:
                j = rng.integers(0, n_seen)
                if j < max_samples:
                    reservoir[j] = item

        for _ in range(n_repeat):
            for node_walk, relation_walk in self.do_walks(length, rng):
                cur_len = 0
                for i in range(len(node_walk)):
                    cur_len = min(cur_len + 1, self._window + 1)
                    if cur_len >= 2:
                        for path_length in range(1, cur_len):
                            rel_seq = tuple(relation_walk[i - path_length:i])
                            path_id = self._path2id[rel_seq]
                            emit((node_walk[i - path_length], node_walk[i], path_id))

        self._id2path = {v: k for k, v in self._path2id.items()}
        print(f"sample(): {n_seen} muestras crudas vistas, {len(reservoir)} conservadas", flush=True)
        return reservoir

    def print_statistics(self):
        print(f"size = {self.node_size}")


def load_a_HIN_from_pandas(edges, print_graph=False):
    """
    edges: DataFrame o lista de DataFrames con columnas:
        source_node, source_class, dest_node, dest_class, edge_class, weight

    El grafo se trata como no dirigido: se agrega automaticamente la
    direccion inversa, con edge_class invertido
    (ej. 'GEN-ass-DIS' -> 'DIS-ass-GEN').
    """

    def reverse(df):
        df = df.rename({"source_node": "dest_node", "dest_node": "source_node",
                         "source_class": "dest_class", "dest_class": "source_class"},
                        axis=1)
        df = df.copy()
        df["edge_class"] = df["edge_class"].map(lambda x: "-".join(reversed(x.split("-"))))
        return df

    print("load graph from edges...", flush=True)
    if isinstance(edges, list):
        edges = pd.concat(edges, sort=False)
    edges = pd.concat([edges, reverse(edges)], sort=False, ignore_index=True)

    g = HIN()
    g.build_from_edges(edges)

    if print_graph:
        g.print_statistics()
    print("finish loading graph!", flush=True)
    return g
