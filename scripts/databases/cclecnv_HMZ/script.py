#!/usr/bin/env python3
"""
cclecnv_HMZ -- script.py
Procesado de CCLE CNV profiles via Harmonizome (gene_attribute_edges.txt.gz).

Logica (fiel a Bioteque):
  1. Leer archivo de edges de Harmonizome
  2. Separar por weight: 1.0 -> CLL-cnu-GEN, -1.0 -> CLL-cnd-GEN
  3. Mapear gene symbol -> ENSG via HGNC
  4. Mapear cell line name -> CVCL via StrippedCellLineName en Model.csv

Columnas del archivo (TSV):
  source | source_desc | source_id | target | target_desc | target_id | weight
  gene_symbol | na | entrez_id | cell_line_name | tissue | -666 | 1.0 or -1.0

Output:
  processed/cclecnv_HMZ/CLL-cnu-GEN.tsv   (n1=CVCL, n2=ENSG) copy number up
  processed/cclecnv_HMZ/CLL-cnd-GEN.tsv                        copy number down
  processed/cclecnv_HMZ/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --symbol2ensg <ruta_symbol2ensg.tsv>
                   --model-csv <ruta_Model.csv>
"""

import argparse
import gzip
import os
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',    required=True)
    parser.add_argument('--out-dir',    required=True)
    parser.add_argument('--symbol2ensg', required=True)
    parser.add_argument('--model-csv',  required=True,
                        help='Ruta a Model.csv de DepMap (tiene StrippedCellLineName y RRID)')
    return parser.parse_args()


def load_symbol2ensg(path: str) -> dict:
    df = pd.read_csv(path, sep='\t', low_memory=False,
                     usecols=['symbol', 'ensembl_gene_id'])
    df = df.dropna(subset=['ensembl_gene_id'])
    mapping = dict(zip(df['symbol'], df['ensembl_gene_id']))
    print(f'[cclecnv_HMZ] Gene symbols mapeados: {len(mapping):,}')
    return mapping


def load_cl2cvcl(model_csv: str) -> dict:
    """Mapea StrippedCellLineName -> CVCL via Model.csv."""
    df = pd.read_csv(model_csv, usecols=['StrippedCellLineName', 'RRID'],
                     low_memory=False)
    df = df.dropna(subset=['RRID'])
    df = df[df['RRID'].str.startswith('CVCL_')]
    mapping = dict(zip(df['StrippedCellLineName'], df['RRID']))
    print(f'[cclecnv_HMZ] Cell lines mapeadas a CVCL: {len(mapping):,}')
    return mapping


def process(raw_dir: str, out_dir: str, symbol2ensg: dict, cl2cvcl: dict):
    os.makedirs(out_dir, exist_ok=True)

    infile = os.path.join(raw_dir, 'gene_attribute_edges.txt.gz')
    print(f'[cclecnv_HMZ] Procesando {infile} ...')

    cnu = set()  # copy number up   (weight == 1)
    cnd = set()  # copy number down (weight == -1)

    skipped_gene   = 0
    skipped_cl     = 0
    skipped_weight = 0
    n_lines        = 0

    with gzip.open(infile, 'rt') as f:
        next(f)  # cabecera 1: column names
        next(f)  # cabecera 2: type descriptions

        for line in f:
            n_lines += 1
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 7:
                continue

            symbol  = cols[0].strip()
            cl_name = cols[3].strip()
            try:
                weight = float(cols[6].strip())
            except ValueError:
                skipped_weight += 1
                continue

            if weight not in (1.0, -1.0):
                skipped_weight += 1
                continue

            ensg = symbol2ensg.get(symbol)
            if ensg is None:
                skipped_gene += 1
                continue

            cvcl = cl2cvcl.get(cl_name)
            if cvcl is None:
                skipped_cl += 1
                continue

            if weight == 1.0:
                cnu.add((cvcl, ensg))
            else:
                cnd.add((cvcl, ensg))

    print(f'[cclecnv_HMZ] Lineas procesadas:        {n_lines:,}')
    print(f'[cclecnv_HMZ] Sin mapeo gene symbol:    {skipped_gene:,}')
    print(f'[cclecnv_HMZ] Sin mapeo cell line:      {skipped_cl:,}')
    print(f'[cclecnv_HMZ] Weight invalido:          {skipped_weight:,}')
    print(f'[cclecnv_HMZ] CLL-cnu-GEN: {len(cnu):,}')
    print(f'[cclecnv_HMZ] CLL-cnd-GEN: {len(cnd):,}')

    # Guardar
    cnu_path = os.path.join(out_dir, 'CLL-cnu-GEN.tsv')
    cnd_path = os.path.join(out_dir, 'CLL-cnd-GEN.tsv')

    with open(cnu_path, 'w') as f:
        f.write('n1\tn2\n')
        for cvcl, ensg in sorted(cnu):
            f.write(f'{cvcl}\t{ensg}\n')

    with open(cnd_path, 'w') as f:
        f.write('n1\tn2\n')
        for cvcl, ensg in sorted(cnd):
            f.write(f'{cvcl}\t{ensg}\n')

    print(f'[cclecnv_HMZ] CLL-cnu-GEN -> {cnu_path}')
    print(f'[cclecnv_HMZ] CLL-cnd-GEN -> {cnd_path}')

    # Stats
    total_db = len(cnu) + len(cnd)
    stats_rows = []
    for rel, s in [('CLL-cnu-GEN', cnu), ('CLL-cnd-GEN', cnd)]:
        stats_rows.append({
            'database':      'cclecnv_HMZ',
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
    print('[cclecnv_HMZ] Procesado completado.')


def main():
    args = parse_args()
    symbol2ensg = load_symbol2ensg(args.symbol2ensg)
    cl2cvcl     = load_cl2cvcl(args.model_csv)
    process(args.raw_dir, args.out_dir, symbol2ensg, cl2cvcl)
    sys.stderr.write('[cclecnv_HMZ] Done!\n')


if __name__ == '__main__':
    main()
