#!/usr/bin/env python3
"""
build_graph_table_v4.py
Cambios respecto a v3:
  - Todas las 17 bases de datos actualizadas a su version mas reciente
    (ver tabla de versiones del proyecto)
  - opentargets: CAMBIO SEMANTICO -- association_overall_direct (v3) ->
    association_overall_indirect (v4). Indirect incluye asociaciones
    propagadas via jerarquia de ontologia EFO/MONDO/HP, alineado con
    la metodologia de Bioteque. Resultado: 283,977 asociaciones GEN-ass-DIS
    (vs el conteo de v3 en association_overall_direct).
  - ctdchemis: se agrego validacion anti-captcha en get_data.sh (no
    afecta el esquema de datos, solo la robustez de la descarga)
  - cclerna_HMZ: fix de bug en script.py (ModelID es columna nombrada
    en el CSV de DepMap 26Q1, no columna de indice sin nombre como
    en versiones anteriores)
  - No se aplico ningun filtro de exclusion de nodos aislados aqui --
    el analisis de componentes conexas se hace DESPUES de construir
    el grafo completo (ver scripts/analysis/), igual que en v3. La
    lista de 8 nodos aislados de v3 no es reutilizable porque el
    cambio direct->indirect en opentargets modifica completamente el
    conjunto de aristas GEN-ass-DIS.
Output en data/processed/integrated/graph_v4/:
  graph_edges.tsv    -- todas las aristas
  graph_summary.tsv  -- resumen por database/relacion
  graph_summary.html -- resumen interactivo + muestra 1000 filas
  graph_stats.tsv    -- metricas globales del grafo
Uso:
  python build_graph_table_v4.py [--base-dir PATH] [--out-dir PATH]
"""
import argparse
import os
import sys
import pandas as pd

BASE_DIR = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity"

RELATIONS = [
    # opentargets -- v4: association_overall_indirect (antes direct en v3)
    ('opentargets',         'GEN-ass-DIS', 'GEN-ass-DIS.tsv', 'n1', 'n2', 'GEN', 'DIS'),
    # hpa_proteome
    ('hpa_proteome',        'GEN-pab-TIS', 'GEN-pab-TIS.tsv', 'n1', 'n2', 'GEN', 'TIS'),
    ('hpa_proteome',        'GEN-pdf-TIS', 'GEN-pdf-TIS.tsv', 'n1', 'n2', 'GEN', 'TIS'),
    # hpa_rna_cons
    ('hpa_rna_cons',        'GEN-upr-TIS', 'GEN-upr-TIS.tsv', 'n1', 'n2', 'GEN', 'TIS'),
    ('hpa_rna_cons',        'GEN-dwr-TIS', 'GEN-dwr-TIS.tsv', 'n1', 'n2', 'GEN', 'TIS'),
    # jensentissuecurated
    ('jensentissuecurated', 'GEN-ass-TIS', 'GEN-ass-TIS.tsv', 'n1', 'n2', 'GEN', 'TIS'),
    # reactome
    ('reactome',            'GEN-ass-PWY', 'GEN-ass-PWY.tsv', 'n1', 'n2', 'GEN', 'PWY'),
    # repohub
    ('repohub',             'CPD-int-GEN', 'CPD-int-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    # coexpressdb
    ('coexpressdb',         'GEN-cex-GEN', 'GEN-cex-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    # omnipath
    ('omnipath',            'GEN-ppi-GEN', 'GEN-ppi-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('omnipath',            'GEN-pho-GEN', 'GEN-pho-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('omnipath',            'GEN-dph-GEN', 'GEN-dph-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    # string
    ('string',              'GEN-ppi-GEN', 'GEN-ppi-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    # dorothea_AB
    ('dorothea_AB',         'GEN-reg-GEN', 'GEN-reg-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_AB',         'GEN-upr-GEN', 'GEN-upr-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_AB',         'GEN-dwr-GEN', 'GEN-dwr-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    # dorothea_CD
    ('dorothea_CD',         'GEN-reg-GEN', 'GEN-reg-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_CD',         'GEN-upr-GEN', 'GEN-upr-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_CD',         'GEN-dwr-GEN', 'GEN-dwr-GEN.tsv', 'n1', 'n2', 'GEN', 'GEN'),
    # cclerna_HMZ
    ('cclerna_HMZ',         'CLL-upr-GEN', 'CLL-upr-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    ('cclerna_HMZ',         'CLL-dwr-GEN', 'CLL-dwr-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    # cclemut_HMZ
    ('cclemut_HMZ',         'CLL-mut-GEN', 'CLL-mut-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    # cclecnv_HMZ
    ('cclecnv_HMZ',         'CLL-cnu-GEN', 'CLL-cnu-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    ('cclecnv_HMZ',         'CLL-cnd-GEN', 'CLL-cnd-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    # achilles_HMZ
    ('achilles_HMZ',        'CLL-bfn-GEN', 'CLL-bfn-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    ('achilles_HMZ',        'CLL-gfn-GEN', 'CLL-gfn-GEN.tsv', 'n1', 'n2', 'CLL', 'GEN'),
    # ctdchemgen (10 relaciones, igual que v3)
    ('ctdchemgen',          'CPD-upe-GEN', 'CPD-upe-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-dwe-GEN', 'CPD-dwe-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-afe-GEN', 'CPD-afe-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-upm-GEN', 'CPD-upm-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-dwm-GEN', 'CPD-dwm-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-afm-GEN', 'CPD-afm-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-upb-GEN', 'CPD-upb-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-dwb-GEN', 'CPD-dwb-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-upa-GEN', 'CPD-upa-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    ('ctdchemgen',          'CPD-dwa-GEN', 'CPD-dwa-GEN.tsv', 'n1', 'n2', 'CPD', 'GEN'),
    # ctdchemis (CPD-DIS)
    ('ctdchemis',           'CPD-cau-DIS', 'CPD-cau-DIS.tsv', 'n1', 'n2', 'CPD', 'DIS'),
    ('ctdchemis',           'CPD-trt-DIS', 'CPD-trt-DIS.tsv', 'n1', 'n2', 'CPD', 'DIS'),
]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-dir', default=BASE_DIR)
    parser.add_argument('--out-dir',  default=None)
    return parser.parse_args()

def load_relation(proc_dir, db, tsv, col1, col2, type1, type2, relation):
    path = os.path.join(proc_dir, db, tsv)
    if not os.path.exists(path):
        print(f"  [WARN] No encontrado: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path, sep='\t', low_memory=False, usecols=[col1, col2])
    df = df.dropna(subset=[col1, col2])
    df = df.rename(columns={col1: 'node1', col2: 'node2'})
    df['node1_type'] = type1
    df['node2_type'] = type2
    df['relation']   = relation
    df['database']   = db
    df['node1'] = df['node1'].astype(str)
    df['node2'] = df['node2'].astype(str)
    return df[['node1', 'node2', 'node1_type', 'node2_type', 'relation', 'database']]

def build_summary(df):
    rows = []
    for (db, rel), group in df.groupby(['database', 'relation']):
        rows.append({
            'Database':     db,
            'Relation':     rel,
            'Node1 Type':   group['node1_type'].iloc[0],
            'Node2 Type':   group['node2_type'].iloc[0],
            'Total Edges':  f"{len(group):,}",
            'Unique Node1': f"{group['node1'].nunique():,}",
            'Unique Node2': f"{group['node2'].nunique():,}",
        })
    summary = pd.DataFrame(rows)
    totals = df.groupby('database').agg(
        Total_Edges=('node1', 'count'),
        Unique_Node1=('node1', 'nunique'),
        Unique_Node2=('node2', 'nunique'),
    ).reset_index()
    return summary, totals

def compute_stats(df):
    all_nodes_n1 = df[['node1', 'node1_type']].rename(columns={'node1': 'node', 'node1_type': 'node_type'})
    all_nodes_n2 = df[['node2', 'node2_type']].rename(columns={'node2': 'node', 'node2_type': 'node_type'})
    all_nodes = pd.concat([all_nodes_n1, all_nodes_n2]).drop_duplicates(subset='node')
    node_counts = all_nodes.groupby('node_type')['node'].count().sort_values(ascending=False)
    stats = {
        'total_edges':          len(df),
        'unique_databases':     df['database'].nunique(),
        'unique_relation_types': df['relation'].nunique(),
        'unique_nodes_total':   len(all_nodes),
        'node_types':           ','.join(node_counts.index.tolist()),
    }
    for ntype, count in node_counts.items():
        stats[f'nodes_{ntype}'] = count
    print("\n=== GRAPH v4 STATS ===")
    print(f"  Total aristas (filas):     {stats['total_edges']:,}")
    print(f"  Databases:                 {stats['unique_databases']}")
    print(f"  Tipos de relacion:         {stats['unique_relation_types']}")
    print(f"  Nodos unicos totales:      {stats['unique_nodes_total']:,}")
    print(f"  Tipos de nodo:             {stats['node_types']}")
    for ntype, count in node_counts.items():
        print(f"    {ntype}: {count:,}")
    print("======================\n")
    return stats, node_counts

def export_stats(stats, node_counts, out_path):
    rows = [{'metric': k, 'value': v} for k, v in stats.items()]
    pd.DataFrame(rows).to_csv(out_path, sep='\t', index=False)
    print(f"[graph] Stats -> {out_path}")

def export_html(df, summary, totals, node_counts, out_path):
    def df_to_dt(data, tid, page_len=20):
        headers  = ''.join(f'<th>{c}</th>' for c in data.columns)
        rows_html = ''.join(
            '<tr>' + ''.join(f'<td>{v}</td>' for v in row) + '</tr>'
            for row in data.values
        )
        return f'''
        <table id="{tid}" class="display compact nowrap" style="width:100%">
          <thead><tr>{headers}</tr></thead>
          <tbody>{rows_html}</tbody>
        </table>
        <script>
        $(document).ready(function() {{
          $('#{tid}').DataTable({{scrollX:true, pageLength:{page_len}}});
        }});
        </script>'''
    sample = df.sample(min(1000, len(df)), random_state=42).sort_values(['database', 'relation'])
    totals_renamed = totals.rename(columns={
        'database': 'Database', 'Total_Edges': 'Total Edges',
        'Unique_Node1': 'Unique Node1', 'Unique_Node2': 'Unique Node2',
    })
    total_edges = len(df)
    total_dbs   = df['database'].nunique()
    total_rels  = df['relation'].nunique()
    total_nodes = (
        pd.concat([
            df[['node1', 'node1_type']].rename(columns={'node1': 'n', 'node1_type': 't'}),
            df[['node2', 'node2_type']].rename(columns={'node2': 'n', 'node2_type': 't'}),
        ]).drop_duplicates(subset='n')['n'].nunique()
    )
    node_badges = ''.join(
        f'<span class="node-{t}">{t}: {c:,}</span> '
        for t, c in node_counts.items()
    )
    html = f'''<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>TFM-NetActivity - Knowledge Graph v4</title>
<link rel="stylesheet"
  href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
<link rel="stylesheet"
  href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
  h1   {{ color: #2c5fa8; }}
  h3   {{ color: #444; margin-top: 30px; }}
  .card {{ background: white; border-radius: 8px; padding: 20px;
           margin-bottom: 20px; box-shadow: 0 2px 6px rgba(0,0,0,0.08); }}
  .badge-count {{ background: #2c5fa8; color: white; border-radius: 12px;
                  padding: 4px 12px; font-size: 1.1em; margin: 0 6px; }}
  .node-GEN {{ background: #d4edda; color: #155724; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
  .node-DIS {{ background: #f8d7da; color: #721c24; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
  .node-TIS {{ background: #fff3cd; color: #856404; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
  .node-CPD {{ background: #cce5ff; color: #004085; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
  .node-PWY {{ background: #e2d9f3; color: #432874; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
  .node-CLL {{ background: #d1ecf1; color: #0c5460; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; margin: 2px; display:inline-block; }}
</style>
</head>
<body>
<div class="card">
  <h1>TFM-NetActivity -- Knowledge Graph v4</h1>
  <p>
    <span class="badge-count">{total_edges:,} aristas totales</span>
    <span class="badge-count">{total_nodes:,} nodos unicos</span>
    <span class="badge-count">{total_dbs} databases</span>
    <span class="badge-count">{total_rels} tipos de relacion</span>
  </p>
  <p>{node_badges}</p>
  <p style="color:#666; font-size:0.85em;">
    Cambios vs v3: las 17 bases actualizadas a su version mas reciente &mdash;
    opentargets cambia de association_overall_direct a association_overall_indirect
  </p>
</div>
<div class="card">
  <h3>Resumen por database y relacion</h3>
  {df_to_dt(summary, 'tbl_summary', 30)}
</div>
<div class="card">
  <h3>Totales por database</h3>
  {df_to_dt(totals_renamed, 'tbl_totals', 20)}
</div>
<div class="card">
  <h3>Muestra de aristas (1,000 filas aleatorias de {total_edges:,})</h3>
  <p style="color:#666; font-size:0.9em;">
    Para el dataset completo usa <code>graph_edges.tsv</code>
  </p>
  {df_to_dt(sample, 'tbl_sample', 25)}
</div>
</body>
</html>'''
    with open(out_path, 'w') as f:
        f.write(html)
    print(f"[graph] HTML -> {out_path}")

def main():
    args     = parse_args()
    base     = args.base_dir
    proc_dir = os.path.join(base, 'data', 'processed', 'databases')
    out_dir  = args.out_dir or os.path.join(base, 'data', 'processed', 'integrated', 'graph_v4')
    os.makedirs(out_dir, exist_ok=True)
    print("=" * 60)
    print(" TFM-NetActivity -- build_graph_table_v4.py")
    print("=" * 60)
    chunks = []
    for (db, rel, tsv, col1, col2, t1, t2) in RELATIONS:
        print(f"[graph] Cargando {db}/{rel}...")
        chunk = load_relation(proc_dir, db, tsv, col1, col2, t1, t2, rel)
        if len(chunk):
            print(f"  -> {len(chunk):,} aristas")
            chunks.append(chunk)
    if not chunks:
        print("[graph] ERROR: No se cargaron datos.")
        sys.exit(1)
    print("[graph] Concatenando...")
    df = pd.concat(chunks, ignore_index=True)
    print(f"[graph] Total aristas: {len(df):,}")
    tsv_path = os.path.join(out_dir, 'graph_edges.tsv')
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"[graph] TSV -> {tsv_path}")
    print(f"[graph] Tamano: {os.path.getsize(tsv_path)/1e6:.1f} MB")
    summary, totals = build_summary(df)
    summary.to_csv(os.path.join(out_dir, 'graph_summary.tsv'), sep='\t', index=False)
    stats, node_counts = compute_stats(df)
    export_stats(stats, node_counts, os.path.join(out_dir, 'graph_stats.tsv'))
    export_html(df, summary, totals, node_counts, os.path.join(out_dir, 'graph_summary.html'))
    print("\n[graph] Completado. Archivos en:", out_dir)
    print(f"  graph_edges.tsv    : {len(df):,} filas")
    print(f"  graph_summary.tsv  : {len(summary)} relaciones")
    print(f"  graph_stats.tsv    : metricas globales")
    print(f"  graph_summary.html : resumen interactivo")

if __name__ == '__main__':
    main()
