#!/usr/bin/env python3
"""
cclemut_HMZ -- script.py
Procesado de CCLE somatic mutations via DepMap.

Logica (fiel a Bioteque):
  1. Leer matrices binarias de mutaciones (damaging y hotspot)
  2. Filtrar IsDefaultEntryForModel == Yes (una entrada por modelo)
  3. Mapear ModelID (ACH-) -> RRID (CVCL_) via Model.csv
  4. Mapear gene symbol -> ENSG via HGNC
  5. Merge de muestras duplicadas (por max)
  6. Merge de genes duplicados (por max)
  7. Extraer pares (CVCL, ENSG) donde valor == 1

Archivos necesarios:
  OmicsSomaticMutationsMatrixDamaging.csv
  OmicsSomaticMutationsMatrixHotspot.csv
  Model.csv -> mapeo ModelID -> RRID

Output:
  processed/cclemut_HMZ/CLL-mut-GEN.tsv   (n1=CVCL, n2=ENSG)
  processed/cclemut_HMZ/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --symbol2ensg <ruta_symbol2ensg.tsv>
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',     required=True)
    parser.add_argument('--out-dir',     required=True)
    parser.add_argument('--symbol2ensg', required=True)
    return parser.parse_args()


def load_symbol2ensg(path: str) -> dict:
    df = pd.read_csv(path, sep='\t', low_memory=False,
                     usecols=['symbol', 'ensembl_gene_id'])
    df = df.dropna(subset=['ensembl_gene_id'])
    mapping = dict(zip(df['symbol'], df['ensembl_gene_id']))
    print(f'[cclemut_HMZ] Gene symbols mapeados: {len(mapping):,}')
    return mapping


def load_model_mapping(raw_dir: str) -> dict:
    path = os.path.join(raw_dir, 'Model.csv')
    df = pd.read_csv(path, usecols=['ModelID', 'RRID'], low_memory=False)
    df = df.dropna(subset=['RRID'])
    df = df[df['RRID'].str.startswith('CVCL_')]
    mapping = dict(zip(df['ModelID'], df['RRID']))
    print(f'[cclemut_HMZ] Modelos mapeados a CVCL: {len(mapping):,}')
    return mapping


def load_mutation_matrix(path: str, model2cvcl: dict, symbol2ensg: dict) -> pd.DataFrame:
    """Carga una matriz binaria de mutaciones y la mapea a CVCL x ENSG."""
    print(f'[cclemut_HMZ] Cargando {os.path.basename(path)} ...')
    df = pd.read_csv(path, index_col=0, low_memory=False)
    print(f'  Shape original: {df.shape}')

    # Filtrar entrada por defecto por modelo
    if 'IsDefaultEntryForModel' in df.columns:
        df = df[df['IsDefaultEntryForModel'] == 'Yes']

    # Usar ModelID como index
    if 'ModelID' in df.columns:
        df.index = df['ModelID']

    # Eliminar columnas de metadata
    meta_cols = ['SequencingID', 'ModelID', 'ModelConditionID',
                 'IsDefaultEntryForModel', 'IsDefaultEntryForMC']
    df = df.drop(columns=[c for c in meta_cols if c in df.columns])

    # Mapear ModelID -> CVCL
    df.index = [model2cvcl.get(x) for x in df.index]
    df = df[df.index.notna()]
    print(f'  Tras mapeo CVCL: {df.shape[0]} lineas celulares')

    # Mapear gene columns: "SYMBOL (ENTREZ)" -> ENSG
    new_cols = []
    for col in df.columns:
        symbol = col.split(' ')[0].strip()
        ensg = symbol2ensg.get(symbol)
        new_cols.append(ensg)
    df.columns = new_cols
    df = df.loc[:, df.columns.notna()]
    print(f'  Tras mapeo ENSG: {df.shape[1]} genes')

    # Merge muestras duplicadas por max (pandas moderno: groupby por index)
    df = df.groupby(df.index).max()

    # Merge genes duplicados por max
    df = df.T.groupby(df.T.index).max().T

    print(f'  Shape final: {df.shape}')
    return df


def process(raw_dir: str, out_dir: str, symbol2ensg: dict, model2cvcl: dict):
    os.makedirs(out_dir, exist_ok=True)

    edges = set()

    for fname in ['OmicsSomaticMutationsMatrixDamaging.csv',
                  'OmicsSomaticMutationsMatrixHotspot.csv']:
        fpath = os.path.join(raw_dir, fname)
        if not os.path.exists(fpath):
            print(f'[cclemut_HMZ] WARN: {fname} no encontrado, saltando.')
            continue

        df = load_mutation_matrix(fpath, model2cvcl, symbol2ensg)

        gns = np.asarray(df.columns)
        for cl, muts in df.iterrows():
            mut_genes = gns[np.array(muts) == 1]
            edges.update(zip([cl] * len(mut_genes), mut_genes))

        print(f'[cclemut_HMZ] Edges acumulados tras {fname}: {len(edges):,}')

    # Guardar
    out_path = os.path.join(out_dir, 'CLL-mut-GEN.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\n')
        for cl, ensg in sorted(edges):
            f.write(f'{cl}\t{ensg}\n')

    print(f'[cclemut_HMZ] CLL-mut-GEN: {len(edges):,} -> {out_path}')

    # Stats
    df_out = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':      'cclemut_HMZ',
        'relation':      'CLL-mut-GEN',
        'gene_id_type':  'ENSG',
        'assoc_id_type': 'CVCL (cell line)',
        'total_assoc':   len(df_out),
        'unique_genes':  df_out['n2'].nunique(),
        'unique_cl':     df_out['n1'].nunique(),
        'db_total':      len(df_out),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[cclemut_HMZ] Procesado completado.')


def main():
    args = parse_args()
    symbol2ensg = load_symbol2ensg(args.symbol2ensg)
    model2cvcl  = load_model_mapping(args.raw_dir)
    process(args.raw_dir, args.out_dir, symbol2ensg, model2cvcl)
    sys.stderr.write('[cclemut_HMZ] Done!\n')


if __name__ == '__main__':
    main()
