#!/usr/bin/env python3
"""
cclerna_HMZ -- script.py
Procesado de CCLE RNA expression via DepMap (OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv).

Logica (fiel a Bioteque / pipeline Harmonizome):
  1. Eliminar genes con expresion cero en mas del 95% de las muestras
  2. Mapear ModelID (ACH-) -> RRID (CVCL_) via Model.csv
  3. Mapear gene symbol -> ENSG via HGNC (symbol2ensg.tsv)
  4. Merge de muestras duplicadas por media
  5. Merge de genes duplicados por media
  6. Quantile normalization por columna (gen)
  7. Escalado por gen con RobustScaler
  8. Top 250 genes up y down por linea celular

Parametros (identicos a Bioteque):
  max_gns = 250

Archivos necesarios:
  OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
  Model.csv
  symbol2ensg.tsv (HGNC)

Output:
  processed/cclerna_HMZ/CLL-upr-GEN.tsv   (n1=CVCL, n2=ENSG)
  processed/cclerna_HMZ/CLL-dwr-GEN.tsv
  processed/cclerna_HMZ/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --symbol2ensg <ruta_symbol2ensg.tsv>
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import RobustScaler
from tqdm import tqdm


MAX_GNS = 250


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',     required=True)
    parser.add_argument('--out-dir',     required=True)
    parser.add_argument('--symbol2ensg', required=True,
                        help='Ruta a data/processed/symbol2ensg.tsv (HGNC)')
    return parser.parse_args()


def load_symbol2ensg(path: str) -> dict:
    """Carga mapeo gene symbol -> ENSG desde HGNC. Solo primer ENSG por symbol."""
    print(f'[cclerna_HMZ] Cargando mapeo symbol->ENSG desde {path} ...')
    df = pd.read_csv(path, sep='\t', low_memory=False,
                     usecols=['symbol', 'ensembl_gene_id'])
    df = df.dropna(subset=['ensembl_gene_id'])
    mapping = dict(zip(df['symbol'], df['ensembl_gene_id']))
    print(f'[cclerna_HMZ] Gene symbols mapeados: {len(mapping):,}')
    return mapping


def load_model_mapping(raw_dir: str) -> dict:
    """Carga mapeo ModelID (ACH-) -> RRID (CVCL_) desde Model.csv."""
    path = os.path.join(raw_dir, 'Model.csv')
    print(f'[cclerna_HMZ] Cargando mapeo ModelID->RRID desde {path} ...')
    df = pd.read_csv(path, usecols=['ModelID', 'RRID'], low_memory=False)
    df = df.dropna(subset=['RRID'])
    df = df[df['RRID'].str.startswith('CVCL_')]
    mapping = dict(zip(df['ModelID'], df['RRID']))
    print(f'[cclerna_HMZ] Modelos mapeados a CVCL: {len(mapping):,}')
    return mapping


def quantile_norm(matrix: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalization por columna (gen). Replica Bioteque."""
    m = np.asarray(matrix)
    rank_mean = np.nanmean(np.sort(m, axis=1), axis=0)
    m_transf = []
    for i in tqdm(range(m.shape[0]), desc='[cclerna_HMZ] Quantile norm'):
        v = rank_mean[rankdata(m[i], 'min') - 1]
        m_transf.append(v)
    return pd.DataFrame(np.asarray(m_transf),
                        index=matrix.index,
                        columns=matrix.columns)


def process(raw_dir: str, out_dir: str, symbol2ensg: dict, model2cvcl: dict):
    os.makedirs(out_dir, exist_ok=True)

    # Cargar expresion
    expr_path = os.path.join(raw_dir, 'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')
    print(f'[cclerna_HMZ] Cargando {expr_path} ...')
    df = pd.read_csv(expr_path, low_memory=False)
    df = df.set_index('ModelID')
    meta_cols = ['SequencingID', 'ModelConditionID', 'IsDefaultEntryForMC', 'IsDefaultEntryForModel']
    df = df.drop(columns=[c for c in meta_cols if c in df.columns])
    print(f'  Shape original: {df.shape}')

    # 1. Eliminar genes con >95% zeros
    df.replace(0, np.nan, inplace=True)
    df.dropna(thresh=int(0.05 * df.shape[0]), axis=1, inplace=True)
    df.replace(np.nan, 0, inplace=True)
    print(f'  Tras filtro 95% zeros: {df.shape}')

    # 2. Mapear ModelID -> CVCL
    df.index = [model2cvcl.get(x) for x in df.index]
    df = df[df.index.notna()]
    print(f'  Tras mapeo CVCL: {df.shape[0]} lineas celulares')

    # 3. Mapear gene columns: "SYMBOL (ENTREZ)" -> ENSG
    new_cols = []
    for col in df.columns:
        symbol = col.split(' ')[0].strip()
        ensg = symbol2ensg.get(symbol)
        new_cols.append(ensg)
    df.columns = new_cols
    df = df.loc[:, df.columns.notna()]
    print(f'  Tras mapeo ENSG: {df.shape[1]} genes')

    # 4. Merge muestras duplicadas por media (pandas moderno)
    df = df.groupby(df.index).mean()

    # 5. Merge genes duplicados por media
    df = df.T.groupby(df.T.index).mean().T
    print(f'  Tras merge duplicados: {df.shape}')

    # 6. Quantile normalization
    df = quantile_norm(df)

    # 7. Escalado por gen con RobustScaler
    print('[cclerna_HMZ] Escalando con RobustScaler ...')
    scaler = RobustScaler()
    df = pd.DataFrame(scaler.fit_transform(df),
                      index=df.index,
                      columns=df.columns)

    # 8. Top 250 up y down por linea celular
    print('[cclerna_HMZ] Calculando upr/dwr por linea celular ...')
    gns = np.asarray(df.columns)
    upr = set()
    dwr = set()

    for cl, expr in df.iterrows():
        expr = np.asarray(expr)
        ixs  = np.argsort(expr)
        upr.update(zip([cl] * MAX_GNS, gns[ixs[-MAX_GNS:]]))
        dwr.update(zip([cl] * MAX_GNS, gns[ixs[:MAX_GNS]]))

    # Guardar
    upr_path = os.path.join(out_dir, 'CLL-upr-GEN.tsv')
    dwr_path = os.path.join(out_dir, 'CLL-dwr-GEN.tsv')

    with open(upr_path, 'w') as f:
        f.write('n1\tn2\n')
        for cl, ensg in sorted(upr):
            f.write(f'{cl}\t{ensg}\n')

    with open(dwr_path, 'w') as f:
        f.write('n1\tn2\n')
        for cl, ensg in sorted(dwr):
            f.write(f'{cl}\t{ensg}\n')

    print(f'[cclerna_HMZ] CLL-upr-GEN: {len(upr):,} -> {upr_path}')
    print(f'[cclerna_HMZ] CLL-dwr-GEN: {len(dwr):,} -> {dwr_path}')

    # Stats
    total_db = len(upr) + len(dwr)
    stats_rows = []
    for rel, s in [('CLL-upr-GEN', upr), ('CLL-dwr-GEN', dwr)]:
        stats_rows.append({
            'database':      'cclerna_HMZ',
            'relation':      rel,
            'gene_id_type':  'ENSG',
            'assoc_id_type': 'CVCL (cell line)',
            'total_assoc':   len(s),
            'unique_genes':  len({x[1] for x in s}),
            'unique_cl':     len({x[0] for x in s}),
            'db_total':      total_db,
        })
    pd.DataFrame(stats_rows).to_csv(
        os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False
    )
    print('[cclerna_HMZ] Procesado completado.')


def main():
    args = parse_args()
    symbol2ensg = load_symbol2ensg(args.symbol2ensg)
    model2cvcl  = load_model_mapping(args.raw_dir)
    process(args.raw_dir, args.out_dir, symbol2ensg, model2cvcl)
    sys.stderr.write('[cclerna_HMZ] Done!\n')


if __name__ == '__main__':
    main()
