#!/usr/bin/env python3
"""
hpa_proteome — script.py
Procesado de Human Protein Atlas IHC data (normal_ihc_data.tsv).

Columnas del archivo fuente:
  Gene | Gene name | Tissue | IHC tissue name | Cell type | Level | Reliability

Lógica:
  1. Filtrar Reliability == 'Uncertain'
  2. pab (protein abundance) = Level == 'High'
  3. pdf (protein deficiency) = Level == 'Low'
  4. Eliminar incongruencias (mismo ENSG en pab y pdf para el mismo tejido)

Mapeo: NO necesario. La columna 'Gene' ya contiene ENSG directamente.

Output:
  processed/hpa_proteome/GEN-pab-TIS.tsv   (n1=ENSG, n2=tissue)
  processed/hpa_proteome/GEN-pdf-TIS.tsv
  processed/hpa_proteome/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
"""

import argparse
import os
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    return parser.parse_args()


def load(raw_dir: str) -> pd.DataFrame:
    path = os.path.join(raw_dir, 'normal_ihc_data.tsv')
    print(f'[hpa_proteome] Cargando {path} ...')
    df = pd.read_csv(path, sep='\t', low_memory=False)
    print(f'  Filas: {len(df):,} | Columnas: {list(df.columns)}')
    return df


def process(df: pd.DataFrame, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    # 1. Filtrar Reliability == 'Uncertain'
    n_before = len(df)
    df = df[df['Reliability'] != 'Uncertain'].copy()
    print(f'[hpa_proteome] Tras filtro Uncertain: {n_before:,} → {len(df):,} filas')

    # 2. Separar por Level
    pab = set(zip(df.loc[df['Level'] == 'High', 'Gene'],
                  df.loc[df['Level'] == 'High', 'Tissue']))

    pdf = set(zip(df.loc[df['Level'] == 'Low', 'Gene'],
                  df.loc[df['Level'] == 'Low', 'Tissue']))

    # 3. Eliminar incongruencias (mismo par ENSG-Tissue en pab y pdf)
    incongruent = pab & pdf
    if incongruent:
        print(f'[hpa_proteome] Incongruencias eliminadas: {len(incongruent):,}')
    pab -= incongruent
    pdf -= incongruent

    # 4. Guardar
    pab_path = os.path.join(out_dir, 'GEN-pab-TIS.tsv')
    pdf_path = os.path.join(out_dir, 'GEN-pdf-TIS.tsv')

    with open(pab_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, tis in sorted(pab):
            f.write(f'{ensg}\t{tis}\n')

    with open(pdf_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, tis in sorted(pdf):
            f.write(f'{ensg}\t{tis}\n')

    print(f'[hpa_proteome] GEN-pab-TIS: {len(pab):,} pares → {pab_path}')
    print(f'[hpa_proteome] GEN-pdf-TIS: {len(pdf):,} pares → {pdf_path}')

    # 5. Stats
    total_db = len(pab) + len(pdf)
    stats_rows = []
    for rel, s in [('GEN-pab-TIS', pab), ('GEN-pdf-TIS', pdf)]:
        stats_rows.append({
            'database':       'hpa_proteome',
            'relation':       rel,
            'gene_id_type':   'ENSG',
            'assoc_id_type':  'Tissue name',
            'total_assoc':    len(s),
            'unique_genes':   len({x[0] for x in s}),
            'unique_targets': len({x[1] for x in s}),
            'db_total':       total_db,
        })

    stats_path = os.path.join(out_dir, 'stats.tsv')
    pd.DataFrame(stats_rows).to_csv(stats_path, sep='\t', index=False)
    print(f'[hpa_proteome] Stats → {stats_path}')


def main():
    args = parse_args()
    df = load(args.raw_dir)
    process(df, args.out_dir)
    sys.stderr.write('[hpa_proteome] Done!\n')


if __name__ == '__main__':
    main()
