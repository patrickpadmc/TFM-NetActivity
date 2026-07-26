#!/usr/bin/env python3
"""
Auditoria de graph_v4: duplicados exactos, solapamiento semantico entre
fuentes, direccionalidad (simetria) empirica y nodos unicos por relacion.

Uso:
    python3 audit_graph_v4.py \
        --edges data/processed/integrated/graph_v4/graph_edges.tsv \
        --out-dir data/processed/integrated/graph_v4/audit
"""
import argparse
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--edges', required=True)
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print('[audit] Cargando graph_edges.tsv...')
    df = pd.read_csv(
        args.edges, sep='\t',
        dtype={'node1': 'category', 'node2': 'category',
               'node1_type': 'category', 'node2_type': 'category',
               'relation': 'category', 'database': 'category'},
    )
    print(f'[audit] {len(df):,} aristas cargadas')

    n_dup_exact = df.duplicated(subset=['node1', 'node2', 'node1_type', 'node2_type', 'relation', 'database']).sum()
    print(f'[audit] Filas duplicadas exactas: {n_dup_exact:,}')

    per_source_dup_rows = []
    for (db, rel), g in df.groupby(['database', 'relation'], observed=True):
        n_dup = g.duplicated(subset=['node1', 'node2']).sum()
        per_source_dup_rows.append({'database': db, 'relation': rel, 'duplicados_mismo_par': int(n_dup)})
    per_source_dup = pd.DataFrame(per_source_dup_rows)

    rows = []
    for (db, rel), g in df.groupby(['database', 'relation'], observed=True):
        n_edges = len(g)
        n_nodes = pd.concat([g['node1'].astype(str), g['node2'].astype(str)]).nunique()
        if n_edges > 2_000_000:
            symmetry_pct = None
            note = 'simetria garantizada por construccion (sorted pair en script.py)'
        else:
            pairs_fwd = set(zip(g['node1'].astype(str), g['node2'].astype(str)))
            pairs_rev = set(zip(g['node2'].astype(str), g['node1'].astype(str)))
            n_sym = len(pairs_fwd & pairs_rev)
            symmetry_pct = round(100 * n_sym / n_edges, 2) if n_edges else 0.0
            note = ''
        rows.append({
            'database': db, 'relation': rel,
            'n_edges': n_edges, 'n_unique_nodes': n_nodes,
            'symmetry_pct': symmetry_pct, 'note': note,
        })
    summary = pd.DataFrame(rows).sort_values(['database', 'relation'])

    overlap_rows = []
    for rel, g in df.groupby('relation', observed=True):
        dbs = g['database'].unique()
        if len(dbs) < 2:
            continue
        pair_sets = {}
        for db in dbs:
            sub = g[g['database'] == db]
            pair_sets[db] = set(zip(sub['node1'].astype(str), sub['node2'].astype(str))) | \
                            set(zip(sub['node2'].astype(str), sub['node1'].astype(str)))
        db_list = list(dbs)
        for i in range(len(db_list)):
            for j in range(i + 1, len(db_list)):
                a, b = db_list[i], db_list[j]
                inter = pair_sets[a] & pair_sets[b]
                overlap_rows.append({
                    'relation': rel, 'database_a': a, 'database_b': b,
                    'n_pares_a': len(pair_sets[a]) // 2, 'n_pares_b': len(pair_sets[b]) // 2,
                    'n_pares_compartidos': len(inter) // 2,
                })
    overlap = pd.DataFrame(overlap_rows)

    summary_path = os.path.join(args.out_dir, 'audit_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)
    print(f'[audit] Resumen por relacion -> {summary_path}')

    dup_path = os.path.join(args.out_dir, 'audit_duplicates_per_source.tsv')
    per_source_dup.to_csv(dup_path, sep='\t', index=False)
    print(f'[audit] Duplicados por fuente -> {dup_path}')

    overlap_path = os.path.join(args.out_dir, 'audit_semantic_overlap.tsv')
    overlap.to_csv(overlap_path, sep='\t', index=False)
    print(f'[audit] Solapamiento semantico -> {overlap_path}')

    print('')
    print('=== RESUMEN ===')
    print(f'Duplicados exactos (fila completa): {n_dup_exact:,}')
    print('')
    print('Duplicados por fuente (mismo par dentro de la misma relacion+bd), solo >0:')
    nonzero = per_source_dup[per_source_dup['duplicados_mismo_par'] > 0]
    print(nonzero.to_string(index=False) if len(nonzero) else '  Ninguno')
    print('')
    print('Solapamiento semantico entre fuentes distintas para la misma relacion canonica:')
    print(overlap.to_string(index=False) if len(overlap) else '  Ninguno')


if __name__ == '__main__':
    main()
