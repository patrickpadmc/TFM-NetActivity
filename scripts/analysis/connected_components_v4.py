#!/usr/bin/env python3
"""
Analisis de componentes conexas de graph_v4.

Calcula el numero de componentes conexas del grafo no dirigido, sin ponderar,
para identificar nodos aislados (o casi aislados) tras la reconstruccion completa
con las 17 bases de datos actualizadas. No se reutiliza el resultado de graph_v3
(8 nodos aislados), ya que el cambio opentargets direct->indirect modifica
por completo el conjunto de aristas GEN-ass-DIS.

Uso:
    python3 connected_components_v4.py \
        --edges data/processed/integrated/graph_v4/graph_edges.tsv \
        --out-dir data/processed/integrated/graph_v4/connectivity
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--edges', required=True)
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print('[conn] Cargando aristas (dtype category para ahorrar memoria)...')
    df = pd.read_csv(
        args.edges,
        sep='\t',
        usecols=['node1', 'node2'],
        dtype={'node1': 'category', 'node2': 'category'},
    )
    print(f'[conn] {len(df):,} aristas cargadas en {time.time()-t0:.1f}s')

    cat1 = df['node1'].cat.categories
    cat2 = df['node2'].cat.categories
    all_ids = pd.Index(cat1).union(cat2)
    id_to_global = {v: i for i, v in enumerate(all_ids)}
    n_nodes = len(all_ids)
    print(f'[conn] Nodos unicos: {n_nodes:,}')

    map1 = np.array([id_to_global[c] for c in cat1])
    map2 = np.array([id_to_global[c] for c in cat2])

    row = map1[df['node1'].cat.codes.to_numpy()]
    col = map2[df['node2'].cat.codes.to_numpy()]
    data = np.ones(len(row), dtype=np.int8)

    print('[conn] Construyendo matriz dispersa...')
    adj = coo_matrix((data, (row, col)), shape=(n_nodes, n_nodes))

    print('[conn] Calculando componentes conexas...')
    n_components, labels = connected_components(csgraph=adj, directed=False)
    print(f'[conn] Componentes conexas: {n_components:,}')

    comp_sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    print('[conn] Distribucion de tamanos (top 10):')
    print(comp_sizes.head(10).to_string())

    isolated_ids = comp_sizes[comp_sizes == 1].index
    small_ids = comp_sizes[comp_sizes < 5].index

    idx_to_id = np.array(all_ids)
    isolated_nodes = idx_to_id[np.isin(labels, isolated_ids)]
    small_mask = np.isin(labels, small_ids)
    small_comp_nodes = idx_to_id[small_mask]
    small_comp_labels = labels[small_mask]

    print(f'[conn] Nodos aislados (componente tamano 1): {len(isolated_nodes):,}')
    print(f'[conn] Nodos en componentes chicas (<5):     {len(small_comp_nodes):,}')

    comp_sizes_path = os.path.join(args.out_dir, 'component_sizes.tsv')
    comp_sizes.rename('n_nodes').rename_axis('component_id').to_csv(comp_sizes_path, sep='\t')
    print(f'[conn] Distribucion de tamanos -> {comp_sizes_path}')

    isolated_path = os.path.join(args.out_dir, 'isolated_nodes.tsv')
    pd.DataFrame({'node_id': isolated_nodes}).to_csv(isolated_path, sep='\t', index=False)
    print(f'[conn] Nodos aislados -> {isolated_path}')

    small_path = os.path.join(args.out_dir, 'small_component_nodes.tsv')
    pd.DataFrame({'node_id': small_comp_nodes, 'component_id': small_comp_labels}).to_csv(
        small_path, sep='\t', index=False
    )
    print(f'[conn] Nodos en componentes chicas -> {small_path}')

    print(f'[conn] Completado en {time.time()-t0:.1f}s')


if __name__ == '__main__':
    main()
