#!/usr/bin/env python3
"""
hpa_rna_cons - script.py
Procesado de Human Protein Atlas RNA tissue consensus (rna_tissue_consensus.tsv).

Logica (fiel a Bioteque):
  1. Construir matriz tejidos x  genes con valores nTPM
  2. Escalar por gen (columna) con RobustScaler
  3. Para cada tejido, tomar top 250 genes mas altos y mas bajos
  4. Filtrar: upr solo si score >= 0.5, dwr solo si score <= -0.5

Parametros (identicos a Bioteque):
  mx_gns      = 250   (maximo genes up/down por tejido)
  min_pos_expr= 0.5   (umbral minimo para up)
  min_neg_expr= -0.5  (umbral maximo para down)

Columnas del archivo fuente:
  Gene | Gene name | Tissue | nTPM

Mapeo: NO necesario. La columna 'Gene' ya contiene ENSG directamente.

Output:
  processed/hpa_rna_cons/GEN-upr-TIS.tsv   (n1=ENSG, n2=tissue)
  processed/hpa_rna_cons/GEN-dwr-TIS.tsv
  processed/hpa_rna_cons/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import RobustScaler


# Parametros (identicos a Bioteque)
MX_GNS       = 250
MIN_POS_EXPR =  0.5
MIN_NEG_EXPR = -0.5


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    return parser.parse_args()


def load(raw_dir: str) -> pd.DataFrame:
    path = os.path.join(raw_dir, 'rna_tissue_consensus.tsv')
    print(f'[hpa_rna_cons] Cargando {path} ...')
    df = pd.read_csv(path, sep='\t', low_memory=False)
    print(f'  Filas: {len(df):,} | Columnas: {list(df.columns)}')
    return df


def process(data: pd.DataFrame, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    # 1. Construir matriz tejidos x genes (ENSG como columnas)
    print('[hpa_rna_cons] Construyendo matriz tejidos x  genes ...')
    matrix = data.pivot_table(
        index='Tissue',
        columns='Gene',
        values='nTPM',
        aggfunc='mean',   # por si hay duplicados
    )
    print(f'  Dimensiones: {matrix.shape[0]} tejidos x {matrix.shape[1]} genes')

    # 2. Escalar por gen (columna) con RobustScaler
    # RobustScaler usa mediana y IQR -> robusto a outliers de expresion
    print('[hpa_rna_cons] Escalando por gen (RobustScaler) ...')
    scaler = RobustScaler()
    scaled = pd.DataFrame(
        scaler.fit_transform(matrix),
        index=matrix.index,
        columns=matrix.columns,
    )

    # 3. Para cada tejido: top 250 up y down con filtro de score
    print('[hpa_rna_cons] Calculando upr/dwr por tejido ...')
    tis_upr = set()
    tis_dwr = set()

    genes = np.asarray(scaled.columns)

    for tissue, row in scaled.iterrows():
        r = np.asarray(row)
        no_null = ~np.isnan(r)

        # Indices ordenados de menor a mayor score (solo genes con valor)
        sorted_ixs = np.argsort(r[no_null])
        valid_ixs  = np.where(no_null)[0]

        # Up: Ultimos MX_GNS genes con score >= MIN_POS_EXPR
        up_ixs = valid_ixs[sorted_ixs[-MX_GNS:]]
        up_ixs = [i for i in up_ixs if r[i] >= MIN_POS_EXPR]

        # Down: primeros MX_GNS genes con score <= MIN_NEG_EXPR
        dw_ixs = valid_ixs[sorted_ixs[:MX_GNS]]
        dw_ixs = [i for i in dw_ixs if r[i] <= MIN_NEG_EXPR]

        tis_upr.update((genes[i], tissue) for i in up_ixs)
        tis_dwr.update((genes[i], tissue) for i in dw_ixs)

    # 4. Guardar
    upr_path = os.path.join(out_dir, 'GEN-upr-TIS.tsv')
    dwr_path = os.path.join(out_dir, 'GEN-dwr-TIS.tsv')

    with open(upr_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, tis in sorted(tis_upr):
            f.write(f'{ensg}\t{tis}\n')

    with open(dwr_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, tis in sorted(tis_dwr):
            f.write(f'{ensg}\t{tis}\n')

    print(f'[hpa_rna_cons] GEN-upr-TIS: {len(tis_upr):,} pares → {upr_path}')
    print(f'[hpa_rna_cons] GEN-dwr-TIS: {len(tis_dwr):,} pares → {dwr_path}')

    # 5. Stats
    total_db = len(tis_upr) + len(tis_dwr)
    stats_rows = []
    for rel, s in [('GEN-upr-TIS', tis_upr), ('GEN-dwr-TIS', tis_dwr)]:
        stats_rows.append({
            'database':       'hpa_rna_cons',
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
    print(f'[hpa_rna_cons] Stats → {stats_path}')


def main():
    args = parse_args()
    data = load(args.raw_dir)
    process(data, args.out_dir)
    sys.stderr.write('[hpa_rna_cons] Done!\n')


if __name__ == '__main__':
    main()
