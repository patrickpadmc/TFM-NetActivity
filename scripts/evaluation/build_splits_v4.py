#!/usr/bin/env python3
"""
Seccion 3 TFM: particiones train/val/test (70/10/20), negativos validos,
grafo de entrenamiento combinado sin fuga para embeddings, y verificacion
de controles anti-fuga, para las 7 tareas primarias de graph_v4.

Decisiones de diseno (confirmadas por Pat, 2026-07-26):
  - Grafo de entrenamiento UNICO combinado (excluye val+test de las 7 tareas
    a la vez), no un grafo separado por tarea.
  - Negativos: 1:1 en train/val, 1:10 en test.
  - Split estricto por nodos (cold-start): solo GEN-ass-DIS.
  - Reglas de preservacion de conectividad, en dos pasos:
      1) cualquier arista positiva donde alguno de sus dos nodos tenga grado
         global 1 (su unica arista en TODO graph_v4) queda forzada a train,
         nunca elegible para val/test.
      2) tras construir el grafo de entrenamiento combinado, cualquier nodo
         que haya quedado con grado 0 (porque TODAS sus aristas, aun con
         grado>=2, cayeron en val/test de una o mas tareas) se "rescata":
         se recupera una de sus aristas retenidas de vuelta a train.

Uso:
    python3 build_splits_v4.py \
        --edges data/processed/integrated/graph_v4/graph_edges.tsv \
        --out-dir data/processed/evaluation/graph_v4
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

SEED = 42
TRAIN_FRAC, VAL_FRAC, TEST_FRAC = 0.70, 0.10, 0.20
NEG_RATIO = {'train': 1, 'val': 1, 'test': 10}

PRIMARY_TASKS = {
    'GEN-ass-DIS': {'relations': ['GEN-ass-DIS'], 'symmetric': False, 't1': 'GEN', 't2': 'DIS'},
    'GEN-ass-TIS': {'relations': ['GEN-ass-TIS'], 'symmetric': False, 't1': 'GEN', 't2': 'TIS'},
    'GEN-ass-PWY': {'relations': ['GEN-ass-PWY'], 'symmetric': False, 't1': 'GEN', 't2': 'PWY'},
    'CPD-trt-DIS': {'relations': ['CPD-trt-DIS'], 'symmetric': False, 't1': 'CPD', 't2': 'DIS'},
    'CPD-int-GEN': {'relations': ['CPD-int-GEN'], 'symmetric': False, 't1': 'CPD', 't2': 'GEN'},
    'GEN-ppi-GEN': {'relations': ['GEN-ppi-GEN'], 'symmetric': True, 't1': 'GEN', 't2': 'GEN'},
    'CLL-mut-GEN': {'relations': ['CLL-mut-GEN'], 'symmetric': False, 't1': 'CLL', 't2': 'GEN'},
}


def log(msg):
    print(f'[splits] {msg}', flush=True)


def canonicalize(df, c1, c2, symmetric):
    if not symmetric:
        return df[c1].copy(), df[c2].copy()
    swap = df[c1] > df[c2]
    a = df[c1].where(~swap, df[c2])
    b = df[c2].where(~swap, df[c1])
    return a, b


def sample_negatives(rng, t1_nodes, t2_nodes, n_needed, forbidden_set, used_set, symmetric, batch_mult=3):
    negs = []
    have = 0
    guard = 0
    while have < n_needed and guard < 300:
        guard += 1
        remaining = n_needed - have
        batch_size = max(int(remaining * batch_mult), 2000)
        a = rng.choice(t1_nodes, size=batch_size)
        b = rng.choice(t2_nodes, size=batch_size)
        if symmetric:
            swap = a > b
            a2 = np.where(swap, b, a)
            b2 = np.where(swap, a, b)
            a, b = a2, b2
            valid = a != b
            a, b = a[valid], b[valid]
        for x, y in zip(a.tolist(), b.tolist()):
            if have >= n_needed:
                break
            pair = (x, y)
            if pair in forbidden_set or pair in used_set:
                continue
            used_set.add(pair)
            negs.append(pair)
            have += 1
    return negs, have >= n_needed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--edges', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    log('Cargando graph_edges.tsv completo...')
    df = pd.read_csv(
        args.edges, sep='\t',
        dtype={'node1': str, 'node2': str, 'node1_type': 'category',
               'node2_type': 'category', 'relation': 'category', 'database': 'category'},
    )
    log(f'{len(df):,} aristas totales cargadas')

    log('Calculando grado global por nodo (todas las relaciones, toda graph_v4)...')
    deg = pd.concat([df['node1'], df['node2']]).value_counts()
    log(f'{len(deg):,} nodos unicos, grado minimo={deg.min()}, nodos con grado global 1: {(deg == 1).sum():,}')

    log('Precomputando universos de nodos por tipo...')
    node_universe = {}
    for t in ['GEN', 'DIS', 'TIS', 'PWY', 'CPD', 'CLL']:
        a = df.loc[df['node1_type'] == t, 'node1']
        b = df.loc[df['node2_type'] == t, 'node2']
        node_universe[t] = np.array(sorted(set(a.unique()) | set(b.unique())))
        log(f'  {t}: {len(node_universe[t]):,} nodos')

    task_positive_sets = {}
    task_negative_sets = {}
    leak_report_lines = []

    for task_idx, (task_name, cfg) in enumerate(PRIMARY_TASKS.items()):
        log(f'=== Tarea: {task_name} ===')
        sub = df[df['relation'].isin(cfg['relations'])].copy()
        n1c, n2c = canonicalize(sub, 'node1', 'node2', cfg['symmetric'])
        pairs = pd.DataFrame({'n1': n1c, 'n2': n2c}).drop_duplicates().reset_index(drop=True)
        n_total = len(pairs)
        log(f'  Pares unicos (tras canonicalizar/deduplicar): {n_total:,} (filas crudas: {len(sub):,})')

        protected_mask = (deg.reindex(pairs['n1']).to_numpy() == 1) | (deg.reindex(pairs['n2']).to_numpy() == 1)
        n_protected = int(protected_mask.sum())
        log(f'  Aristas protegidas (grado global 1 en algun extremo, forzadas a train): {n_protected}')

        protected = pairs[protected_mask]
        free = pairs[~protected_mask].sample(frac=1.0, random_state=SEED).reset_index(drop=True)

        n_val_target = round(n_total * VAL_FRAC)
        n_test_target = round(n_total * TEST_FRAC)
        n_val = min(n_val_target, len(free))
        n_test = min(n_test_target, len(free) - n_val)

        val = free.iloc[:n_val].reset_index(drop=True)
        test = free.iloc[n_val:n_val + n_test].reset_index(drop=True)
        train_free = free.iloc[n_val + n_test:]
        train = pd.concat([protected, train_free], ignore_index=True)

        log(f'  Split inicial -> train={len(train):,} ({len(train)/n_total:.1%}) '
            f'val={len(val):,} ({len(val)/n_total:.1%}) test={len(test):,} ({len(test)/n_total:.1%})')

        pos_all = set(zip(pairs['n1'], pairs['n2']))
        task_positive_sets[task_name] = {
            'train': set(zip(train['n1'], train['n2'])),
            'val': set(zip(val['n1'], val['n2'])),
            'test': set(zip(test['n1'], test['n2'])),
            'all': pos_all,
            'symmetric': cfg['symmetric'],
            'protected_pairs': set(zip(protected['n1'], protected['n2'])),
        }

        t1_nodes = node_universe[cfg['t1']]
        t2_nodes = node_universe[cfg['t2']]
        log(f'  Pool de negativos: {len(t1_nodes):,} x {len(t2_nodes):,}')

        used_negatives = set()
        neg_splits = {}
        rng_task = np.random.default_rng(SEED + task_idx)
        for split_name, split_set in [('train', task_positive_sets[task_name]['train']),
                                        ('val', task_positive_sets[task_name]['val']),
                                        ('test', task_positive_sets[task_name]['test'])]:
            n_needed = len(split_set) * NEG_RATIO[split_name]
            negs, ok = sample_negatives(rng_task, t1_nodes, t2_nodes, n_needed, pos_all, used_negatives, cfg['symmetric'])
            if not ok:
                log(f'  AVISO: solo se generaron {len(negs)}/{n_needed} negativos para {task_name}/{split_name}')
            neg_splits[split_name] = negs
            log(f'  Negativos {split_name}: {len(negs):,} (ratio objetivo 1:{NEG_RATIO[split_name]})')

        task_negative_sets[task_name] = neg_splits

    log('=== Split estricto (cold-start por nodos) para GEN-ass-DIS ===')
    strict_task = 'GEN-ass-DIS'
    cfg = PRIMARY_TASKS[strict_task]
    sub = df[df['relation'].isin(cfg['relations'])].copy()
    pairs = sub[['node1', 'node2']].drop_duplicates().rename(columns={'node1': 'n1', 'node2': 'n2'}).reset_index(drop=True)
    genes = pairs['n1'].unique()
    gene_deg_global = deg.reindex(genes)
    protected_genes = set(gene_deg_global[gene_deg_global == 1].index)
    free_genes = np.array(sorted(set(genes) - protected_genes))
    rng_g = np.random.default_rng(SEED)
    rng_g.shuffle(free_genes)
    n_genes = len(genes)
    n_val_g = min(round(n_genes * VAL_FRAC), len(free_genes))
    n_test_g = min(round(n_genes * TEST_FRAC), len(free_genes) - n_val_g)
    val_genes = set(free_genes[:n_val_g])
    test_genes = set(free_genes[n_val_g:n_val_g + n_test_g])
    train_genes = set(genes) - val_genes - test_genes

    train_s = pairs[pairs['n1'].isin(train_genes)].reset_index(drop=True)
    val_s = pairs[pairs['n1'].isin(val_genes)].reset_index(drop=True)
    test_s = pairs[pairs['n1'].isin(test_genes)].reset_index(drop=True)
    log(f'  Genes: train={len(train_genes):,} val={len(val_genes):,} test={len(test_genes):,} (protegidos={len(protected_genes)})')
    log(f'  Aristas resultantes: train={len(train_s):,} val={len(val_s):,} test={len(test_s):,}')

    strict_dir = os.path.join(args.out_dir, strict_task)
    os.makedirs(strict_dir, exist_ok=True)
    train_s.to_csv(os.path.join(strict_dir, 'train_strict.tsv'), sep='\t', index=False)
    val_s.to_csv(os.path.join(strict_dir, 'val_strict.tsv'), sep='\t', index=False)
    test_s.to_csv(os.path.join(strict_dir, 'test_strict.tsv'), sep='\t', index=False)

    dis_nodes = node_universe['DIS']
    pos_all_strict = set(zip(pairs['n1'], pairs['n2']))
    rng_gs = np.random.default_rng(SEED + 100)
    used_neg_strict = set()
    strict_neg_counts = {}
    for split_name, split_df, genes_split in [('train', train_s, train_genes), ('val', val_s, val_genes), ('test', test_s, test_genes)]:
        n_needed = len(split_df) * NEG_RATIO[split_name]
        genes_arr = np.array(sorted(genes_split)) if genes_split else np.array([])
        if len(genes_arr) == 0 or n_needed == 0:
            negs = []
        else:
            negs, ok = sample_negatives(rng_gs, genes_arr, dis_nodes, n_needed, pos_all_strict, used_neg_strict, False)
            if not ok:
                log(f'  AVISO: solo se generaron {len(negs)}/{n_needed} negativos estrictos para {split_name}')
        strict_neg_counts[split_name] = len(negs)
        neg_df = pd.DataFrame(negs, columns=['n1', 'n2'])
        neg_df.to_csv(os.path.join(strict_dir, f'negatives_{split_name}_strict.tsv'), sep='\t', index=False)
        log(f'  Negativos estrictos {split_name}: {len(negs):,} (objetivo {n_needed})')

    def build_exclude_mask():
        forbidden = {}
        for tname, tcfg in PRIMARY_TASKS.items():
            s = task_positive_sets[tname]
            forb = s['val'] | s['test']
            for rel in tcfg['relations']:
                forbidden.setdefault(rel, set()).update(forb)
        mask = pd.Series(False, index=df.index)
        for rel, forb_set in forbidden.items():
            rel_mask = df['relation'] == rel
            rel_rows = df.loc[rel_mask, ['node1', 'node2']]
            pairs_fwd = list(zip(rel_rows['node1'], rel_rows['node2']))
            pairs_rev = list(zip(rel_rows['node2'], rel_rows['node1']))
            hit = np.array([p in forb_set or r in forb_set for p, r in zip(pairs_fwd, pairs_rev)])
            mask.loc[rel_rows.index[hit]] = True
        return mask

    log('=== Construyendo grafo de entrenamiento combinado (sin val/test de las 7 tareas) ===')
    exclude_mask = build_exclude_mask()
    train_graph = df.loc[~exclude_mask]
    log(f'  Aristas excluidas (primer intento): {int(exclude_mask.sum()):,}')
    log(f'  Grafo de entrenamiento combinado (primer intento): {len(train_graph):,} aristas')

    full_node_set = set(deg.index)
    graph_node_set = set(train_graph['node1'].unique()) | set(train_graph['node2'].unique())
    stranded = full_node_set - graph_node_set
    log(f'  Nodos que quedarian con grado 0 en el grafo de embeddings: {len(stranded):,}')

    n_rescued = 0
    if stranded:
        log('  Aplicando paso de rescate: recuperando una arista retenida por cada nodo aislado...')
        for node in stranded:
            rescued_this_node = False
            for tname in PRIMARY_TASKS:
                if rescued_this_node:
                    break
                s = task_positive_sets[tname]
                for split_name in ('val', 'test'):
                    candidate = None
                    for pair in s[split_name]:
                        if pair[0] == node or pair[1] == node:
                            candidate = pair
                            break
                    if candidate is not None:
                        s[split_name].discard(candidate)
                        s['train'].add(candidate)
                        n_rescued += 1
                        rescued_this_node = True
                        break
        log(f'  Aristas rescatadas: {n_rescued}')

        exclude_mask = build_exclude_mask()
        train_graph = df.loc[~exclude_mask]
        log(f'  Aristas excluidas (tras rescate): {int(exclude_mask.sum()):,}')
        log(f'  Grafo de entrenamiento combinado (tras rescate): {len(train_graph):,} aristas')

        graph_node_set = set(train_graph['node1'].unique()) | set(train_graph['node2'].unique())
        stranded_after = full_node_set - graph_node_set
        log(f'  Nodos que siguen con grado 0 tras el rescate: {len(stranded_after):,}')
        if stranded_after:
            log(f'  AVISO: {len(stranded_after)} nodos no se pudieron rescatar')
    else:
        log('  No hay nodos en riesgo; no se requiere rescate.')

    n_excluded = int(exclude_mask.sum())

    log('Guardando archivos finales de particiones (reflejando el rescate)...')
    summary_rows = []
    for task_name, cfg in PRIMARY_TASKS.items():
        s = task_positive_sets[task_name]
        task_dir = os.path.join(args.out_dir, task_name)
        os.makedirs(task_dir, exist_ok=True)
        pd.DataFrame(sorted(s['train']), columns=['n1', 'n2']).to_csv(os.path.join(task_dir, 'train.tsv'), sep='\t', index=False)
        pd.DataFrame(sorted(s['val']), columns=['n1', 'n2']).to_csv(os.path.join(task_dir, 'val.tsv'), sep='\t', index=False)
        pd.DataFrame(sorted(s['test']), columns=['n1', 'n2']).to_csv(os.path.join(task_dir, 'test.tsv'), sep='\t', index=False)

        neg = task_negative_sets[task_name]
        for split_name in ('train', 'val', 'test'):
            neg_df = pd.DataFrame(neg[split_name], columns=['n1', 'n2'])
            neg_df.to_csv(os.path.join(task_dir, f'negatives_{split_name}.tsv'), sep='\t', index=False)

        summary_rows.append({
            'task': task_name, 'n_total_pos': len(s['all']), 'n_protected': len(s['protected_pairs']),
            'n_train': len(s['train']), 'n_val': len(s['val']), 'n_test': len(s['test']),
            'n_neg_train': len(neg['train']), 'n_neg_val': len(neg['val']), 'n_neg_test': len(neg['test']),
        })

    train_graph_path = os.path.join(args.out_dir, 'graph_edges_train_only.tsv')
    train_graph.to_csv(train_graph_path, sep='\t', index=False)
    log(f'  Guardado -> {train_graph_path}')

    log('=== Verificando conectividad del grafo de entrenamiento combinado (final) ===')
    cat1 = train_graph['node1'].astype('category')
    cat2 = train_graph['node2'].astype('category')
    all_ids = pd.Index(cat1.cat.categories).union(cat2.cat.categories)
    id_to_idx = {v: i for i, v in enumerate(all_ids)}
    row_idx = train_graph['node1'].map(id_to_idx).to_numpy()
    col_idx = train_graph['node2'].map(id_to_idx).to_numpy()
    n_nodes_g = len(all_ids)
    adj = coo_matrix((np.ones(len(row_idx), dtype=np.int8), (row_idx, col_idx)), shape=(n_nodes_g, n_nodes_g))
    n_comp, labels = connected_components(csgraph=adj, directed=False)
    comp_sizes = pd.Series(labels).value_counts()
    largest = int(comp_sizes.max())
    n_singleton_components = int((comp_sizes == 1).sum())
    log(f'  Nodos totales en grafo de entrenamiento: {n_nodes_g:,} (de {len(full_node_set):,} en graph_v4)')
    log(f'  Componentes conexas: {n_comp:,}; componente principal: {largest:,} nodos ({largest/n_nodes_g:.4%})')
    log(f'  Nodos completamente aislados (componente de tamano 1): {n_singleton_components:,}')

    log('=== Verificacion anti-fuga (final, tras rescate) ===')
    all_ok = True
    for task_name3, s in task_positive_sets.items():
        tr, va, te = s['train'], s['val'], s['test']
        overlap_tv = tr & va
        overlap_tt = tr & te
        overlap_vt = va & te
        ok1 = len(overlap_tv) == 0 and len(overlap_tt) == 0 and len(overlap_vt) == 0
        all_ok = all_ok and ok1
        line = f'{task_name3}: solape positivos train/val={len(overlap_tv)} train/test={len(overlap_tt)} val/test={len(overlap_vt)} -> {"OK" if ok1 else "FALLO"}'
        log('  ' + line)
        leak_report_lines.append(line)

        neg = task_negative_sets[task_name3]
        neg_all = set(neg['train']) | set(neg['val']) | set(neg['test'])
        collision = neg_all & s['all']
        ok2 = len(collision) == 0
        all_ok = all_ok and ok2
        line2 = f'{task_name3}: negativos que coinciden con positivos reales={len(collision)} -> {"OK" if ok2 else "FALLO"}'
        log('  ' + line2)
        leak_report_lines.append(line2)

        neg_overlap_tv = set(neg['train']) & set(neg['val'])
        neg_overlap_tt = set(neg['train']) & set(neg['test'])
        neg_overlap_vt = set(neg['val']) & set(neg['test'])
        ok3 = len(neg_overlap_tv) == 0 and len(neg_overlap_tt) == 0 and len(neg_overlap_vt) == 0
        all_ok = all_ok and ok3
        line3 = f'{task_name3}: negativos duplicados entre splits train/val={len(neg_overlap_tv)} train/test={len(neg_overlap_tt)} val/test={len(neg_overlap_vt)} -> {"OK" if ok3 else "FALLO"}'
        log('  ' + line3)
        leak_report_lines.append(line3)

        prot = s['protected_pairs']
        prot_in_val_test = prot & (va | te)
        ok4 = len(prot_in_val_test) == 0
        all_ok = all_ok and ok4
        line4 = f'{task_name3}: aristas protegidas que terminaron en val/test={len(prot_in_val_test)} -> {"OK" if ok4 else "FALLO"}'
        log('  ' + line4)
        leak_report_lines.append(line4)

    ok5 = n_singleton_components == 0 and (n_nodes_g == len(full_node_set))
    all_ok = all_ok and ok5
    line5 = f'grafo_embeddings: nodos con grado 0 = {len(full_node_set) - n_nodes_g:,}, componentes singleton = {n_singleton_components} -> {"OK" if ok5 else "FALLO"}'
    log('  ' + line5)
    leak_report_lines.append(line5)

    log(f'=== VERIFICACION GLOBAL: {"TODO OK" if all_ok else "HAY FALLOS - REVISAR"} ===')

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.out_dir, 'splits_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    log(f'Resumen -> {summary_path}')

    manifest = {
        'seed': SEED,
        'split_ratios': {'train': TRAIN_FRAC, 'val': VAL_FRAC, 'test': TEST_FRAC},
        'negative_ratios': NEG_RATIO,
        'connectivity_rules': [
            'aristas con algun extremo de grado global 1 forzadas a train (paso 1)',
            'nodos que quedan con grado 0 en el grafo combinado se rescatan recuperando una arista retenida (paso 2)',
        ],
        'strict_split_tasks': ['GEN-ass-DIS'],
        'combined_embedding_graph': True,
        'all_checks_passed': bool(all_ok),
        'n_edges_train_graph': int(len(train_graph)),
        'n_edges_excluded_from_train_graph': n_excluded,
        'n_edges_rescued': n_rescued,
        'connected_components_train_graph': int(n_comp),
        'largest_component_train_graph': largest,
        'singleton_nodes_train_graph': n_singleton_components,
        'strict_split_negative_counts': strict_neg_counts,
    }
    with open(os.path.join(args.out_dir, 'manifest_splits.json'), 'w') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    with open(os.path.join(args.out_dir, 'leakage_report.txt'), 'w') as f:
        f.write('\n'.join(leak_report_lines) + '\n')
        f.write(f'\nVERIFICACION GLOBAL: {"TODO OK" if all_ok else "HAY FALLOS"}\n')

    log('Completado.')


if __name__ == '__main__':
    main()
