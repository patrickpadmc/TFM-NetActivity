#!/usr/bin/env python3
"""
build_graph_table.py
Construye la tabla base para el grafo de conocimiento biologico.

Output:
  data/processed/graph/graph_edges.tsv     -- todas las aristas
  data/processed/graph/graph_summary.html  -- resumen + muestra 1000 filas

Estructura de graph_edges.tsv:
  node1 | node2 | node1_type | node2_type | relation | database

Solo se incluyen relaciones donde al menos un nodo es GEN.

Tipos de nodo:
  GEN -- gen (ENSG)
  DIS -- enfermedad (MESH, EFO, DOID, ORPHA, HP, OMIM, UMLS)
  TIS -- tejido (nombre libre o BTO)
  CPD -- compuesto (CTD ID o pert_iname)
  PWY -- pathway (R-HSA-)
  CLL -- linea celular (CVCL)

Uso:
  python build_graph_table.py [--base-dir PATH] [--out-dir PATH]
"""

import argparse
import os
import sys
import glob
import pandas as pd
import numpy as np


BASE_DIR = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity"

# Definicion de cada relacion:
# (database, relation, file_n1, file_n2, node1_type, node2_type)
# file_n1/n2: nombre de columna en el TSV procesado
RELATIONS = [
    # opentargets
    ('opentargets',         'GEN-ass-DIS', 'GEN-ass-DIS.tsv',     'n1', 'n2', 'GEN', 'DIS'),
    # hpa_proteome
    ('hpa_proteome',        'GEN-pab-TIS', 'GEN-pab-TIS.tsv',     'n1', 'n2', 'GEN', 'TIS'),
    ('hpa_proteome',        'GEN-pdf-TIS', 'GEN-pdf-TIS.tsv',     'n1', 'n2', 'GEN', 'TIS'),
    # hpa_rna_cons
    ('hpa_rna_cons',        'GEN-upr-TIS', 'GEN-upr-TIS.tsv',     'n1', 'n2', 'GEN', 'TIS'),
    ('hpa_rna_cons',        'GEN-dwr-TIS', 'GEN-dwr-TIS.tsv',     'n1', 'n2', 'GEN', 'TIS'),
    # ctddisease
    ('ctddisease',          'GEN-ass-DIS', 'GEN-ass-DIS.tsv',     'n1', 'n2', 'GEN', 'DIS'),
    # jensentissuecurated
    ('jensentissuecurated', 'GEN-ass-TIS', 'GEN-ass-TIS.tsv',     'n1', 'n2', 'GEN', 'TIS'),
    # reactome
    ('reactome',            'GEN-ass-PWY', 'GEN-ass-PWY.tsv',     'n1', 'n2', 'GEN', 'PWY'),
    # repohub
    ('repohub',             'CPD-int-GEN', 'CPD-int-GEN.tsv',     'n1', 'n2', 'CPD', 'GEN'),
    # coexpressdb
    ('coexpressdb',         'GEN-cex-GEN', 'GEN-cex-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    # omnipath
    ('omnipath',            'GEN-ppi-GEN', 'GEN-ppi-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('omnipath',            'GEN-pho-GEN', 'GEN-pho-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('omnipath',            'GEN-dph-GEN', 'GEN-dph-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    # string
    ('string',              'GEN-ppi-GEN', 'GEN-ppi-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    # dorothea_AB
    ('dorothea_AB',         'GEN-reg-GEN', 'GEN-reg-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_AB',         'GEN-upr-GEN', 'GEN-upr-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_AB',         'GEN-dwr-GEN', 'GEN-dwr-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    # dorothea_CD
    ('dorothea_CD',         'GEN-reg-GEN', 'GEN-reg-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_CD',         'GEN-upr-GEN', 'GEN-upr-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    ('dorothea_CD',         'GEN-dwr-GEN', 'GEN-dwr-GEN.tsv',     'n1', 'n2', 'GEN', 'GEN'),
    # cclerna_HMZ
    ('cclerna_HMZ',         'CLL-upr-GEN', 'CLL-upr-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    ('cclerna_HMZ',         'CLL-dwr-GEN', 'CLL-dwr-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    # cclemut_HMZ
    ('cclemut_HMZ',         'CLL-mut-GEN', 'CLL-mut-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    # cclecnv_HMZ
    ('cclecnv_HMZ',         'CLL-cnu-GEN', 'CLL-cnu-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    ('cclecnv_HMZ',         'CLL-cnd-GEN', 'CLL-cnd-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    # achilles_HMZ
    ('achilles_HMZ',        'CLL-bfn-GEN', 'CLL-bfn-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
    ('achilles_HMZ',        'CLL-gfn-GEN', 'CLL-gfn-GEN.tsv',     'n1', 'n2', 'CLL', 'GEN'),
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-dir', default=BASE_DIR)
    parser.add_argument('--out-dir',  default=None)
    return parser.parse_args()


def load_relation(proc_dir: str, db: str, tsv: str,
                  col1: str, col2: str,
                  type1: str, type2: str,
                  relation: str) -> pd.DataFrame:
    """Carga un archivo procesado y lo formatea como aristas del grafo."""
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

    # Castear a string para consistencia
    df['node1'] = df['node1'].astype(str)
    df['node2'] = df['node2'].astype(str)

    return df[['node1', 'node2', 'node1_type', 'node2_type', 'relation', 'database']]


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Genera tabla resumen por database y relacion."""
    rows = []
    for (db, rel), group in df.groupby(['database', 'relation']):
        n1_type = group['node1_type'].iloc[0]
        n2_type = group['node2_type'].iloc[0]
        rows.append({
            'Database':     db,
            'Relation':     rel,
            'Node1 Type':   n1_type,
            'Node2 Type':   n2_type,
            'Total Edges':  f"{len(group):,}",
            'Unique Node1': f"{group['node1'].nunique():,}",
            'Unique Node2': f"{group['node2'].nunique():,}",
        })
    summary = pd.DataFrame(rows)

    # Totales por database
    totals = df.groupby('database').agg(
        Total_Edges=('node1', 'count'),
        Unique_Node1=('node1', 'nunique'),
        Unique_Node2=('node2', 'nunique'),
    ).reset_index()

    return summary, totals


def export_html(df: pd.DataFrame, summary: pd.DataFrame,
                totals: pd.DataFrame, out_path: str):
    """Genera HTML con resumen, totales y muestra de 1000 filas."""

    def df_to_dt(data: pd.DataFrame, tid: str, page_len: int = 20) -> str:
        headers = ''.join(f'<th>{c}</th>' for c in data.columns)
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

    # Muestra de 1000 filas
    sample = df.sample(min(1000, len(df)), random_state=42).sort_values(
        ['database', 'relation']
    )

    # Renombrar totals fuera del f-string para evitar problema con {{ en expresiones
    totals_renamed = totals.rename(columns={
        'database':     'Database',
        'Total_Edges':  'Total Edges',
        'Unique_Node1': 'Unique Node1',
        'Unique_Node2': 'Unique Node2',
    })

    # Conteos globales
    total_edges  = len(df)
    total_genes  = df[df['node1_type']=='GEN']['node1'].nunique()
    total_genes += df[df['node2_type']=='GEN']['node2'].nunique()
    total_dbs    = df['database'].nunique()
    total_rels   = df['relation'].nunique()

    html = f'''<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>TFM-NetActivity - Knowledge Graph Edges</title>
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
               border-radius: 4px; font-size: 0.85em; }}
  .node-DIS {{ background: #f8d7da; color: #721c24; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; }}
  .node-TIS {{ background: #fff3cd; color: #856404; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; }}
  .node-CPD {{ background: #cce5ff; color: #004085; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; }}
  .node-PWY {{ background: #e2d9f3; color: #432874; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; }}
  .node-CLL {{ background: #d1ecf1; color: #0c5460; padding: 2px 8px;
               border-radius: 4px; font-size: 0.85em; }}
</style>
</head>
<body>
<div class="card">
  <h1>TFM-NetActivity -- Knowledge Graph Edges</h1>
  <p>
    <span class="badge-count">{total_edges:,} aristas totales</span>
    <span class="badge-count">{total_dbs} databases</span>
    <span class="badge-count">{total_rels} relaciones</span>
  </p>
  <p>
    <span class="node-GEN">GEN</span>
    <span class="node-DIS">DIS</span>
    <span class="node-TIS">TIS</span>
    <span class="node-CPD">CPD</span>
    <span class="node-PWY">PWY</span>
    <span class="node-CLL">CLL</span>
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
    proc_dir = f"{base}/data/processed"
    out_dir  = args.out_dir or f"{proc_dir}/graph"
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 60)
    print(" TFM-NetActivity -- build_graph_table.py")
    print("=" * 60)

    chunks = []
    for (db, rel, tsv, col1, col2, t1, t2) in RELATIONS:
        print(f"[graph] Cargando {db}/{rel}...")
        chunk = load_relation(proc_dir, db, tsv, col1, col2, t1, t2, rel)
        if len(chunk) > 0:
            print(f"  -> {len(chunk):,} aristas")
            chunks.append(chunk)

    if not chunks:
        print("[graph] ERROR: No se cargaron datos.")
        sys.exit(1)

    print("[graph] Concatenando...")
    df = pd.concat(chunks, ignore_index=True)
    print(f"[graph] Total aristas: {len(df):,}")

    # Guardar TSV
    tsv_path = os.path.join(out_dir, 'graph_edges.tsv')
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"[graph] TSV -> {tsv_path}")
    print(f"[graph] Tamano: {os.path.getsize(tsv_path)/1e6:.1f} MB")

    # Resumen y HTML
    summary, totals = build_summary(df)
    summary_path = os.path.join(out_dir, 'graph_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)
    print(f"[graph] Summary TSV -> {summary_path}")

    html_path = os.path.join(out_dir, 'graph_summary.html')
    export_html(df, summary, totals, html_path)

    print("\n[graph] Completado. Archivos en:", out_dir)
    print(f"  graph_edges.tsv   : {len(df):,} filas")
    print(f"  graph_summary.tsv : {len(summary)} relaciones")
    print(f"  graph_summary.html: resumen + muestra 1000 filas")


if __name__ == '__main__':
    main()
