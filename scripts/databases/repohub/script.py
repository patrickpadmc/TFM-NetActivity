#!/usr/bin/env python3
"""
repohub -- script.py
Procesado de Drug Repurposing Hub (repo-drug-annotation-20200324.txt).

Logica (fiel a Bioteque):
  1. Leer archivo saltando lineas con !
  2. Mapear gene symbol -> ENSG via HGNC complete set
  3. Construir pares (compound, ENSG)
  4. Eliminar outliers por zscore > 3 en grado de conectividad

Columnas del archivo:
  pert_iname | clinical_phase | moa | target | disease_area | indication

Mapeo: gene symbol -> ENSG via data/processed/symbol2ensg.tsv
  (HGNC complete set, columnas: symbol, ensembl_gene_id)

Output:
  processed/repohub/CPD-int-GEN.tsv   (n1=compound, n2=ENSG)
  processed/repohub/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --symbol2ensg <ruta_symbol2ensg.tsv>
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import zscore


ZSCORE_CUTOFF = 3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',     required=True)
    parser.add_argument('--out-dir',     required=True)
    parser.add_argument('--symbol2ensg', required=True,
                        help='Ruta a data/processed/symbol2ensg.tsv (HGNC)')
    return parser.parse_args()


def load_symbol2ensg(path: str) -> dict:
    """Carga mapeo gene symbol -> ENSG desde HGNC complete set."""
    print(f'[repohub] Cargando mapeo gene symbol->ENSG desde {path} ...')
    df = pd.read_csv(path, sep='\t', low_memory=False,
                     usecols=['symbol', 'prev_symbol', 'alias_symbol', 'ensembl_gene_id'])
    df = df.dropna(subset=['ensembl_gene_id'])

    mapping = {}

    # Simbolo principal
    for _, row in df.iterrows():
        ensg = row['ensembl_gene_id'].strip()
        sym  = row['symbol'].strip()
        if sym:
            if sym not in mapping:
                mapping[sym] = set()
            mapping[sym].add(ensg)

        # Simbolos previos (gene renaming)
        prev = str(row.get('prev_symbol', '') or '')
        for s in prev.split('|'):
            s = s.strip()
            if s and s not in mapping:
                mapping[s] = set()
            if s:
                mapping[s].add(ensg)

        # Alias
        alias = str(row.get('alias_symbol', '') or '')
        for s in alias.split('|'):
            s = s.strip()
            if s and s not in mapping:
                mapping[s] = set()
            if s:
                mapping[s].add(ensg)

    print(f'[repohub] Gene symbols mapeados: {len(mapping):,}')
    return mapping


def detect_outliers_zscore(nodes: list, values: list, cutoff: float = ZSCORE_CUTOFF) -> set:
    """Detecta outliers por zscore > cutoff. Replica detect_outliers de Bioteque."""
    if len(values) < 3:
        return set()
    scores = zscore(values)
    return {n for n, s in zip(nodes, scores) if s > cutoff}


def find_outliers_from_edges(edges: set) -> tuple:
    """Replica find_outliers_from_edges de Bioteque (directed, zscore)."""
    n1_counts = Counter(e[0] for e in edges)
    n2_counts = Counter(e[1] for e in edges)

    n1_outliers = detect_outliers_zscore(
        list(n1_counts.keys()), list(n1_counts.values())
    )
    n2_outliers = detect_outliers_zscore(
        list(n2_counts.keys()), list(n2_counts.values())
    )
    return n1_outliers, n2_outliers


def load_data(raw_dir: str) -> pd.DataFrame:
    infile = os.path.join(raw_dir, 'repo-drug-annotation-20200324.txt')
    print(f'[repohub] Cargando {infile} ...')
    rows = []
    header = None
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue
            cols = line.rstrip('\n').split('\t')
            if header is None:
                header = cols
                continue
            rows.append(cols)
    df = pd.DataFrame(rows, columns=header)
    print(f'[repohub] Filas: {len(df):,} | Columnas: {list(df.columns)}')
    return df


def process(df: pd.DataFrame, out_dir: str, symbol2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)

    edges = set()
    skipped_no_target  = 0
    skipped_no_mapping = 0

    for _, row in df.iterrows():
        compound = row.get('pert_iname', '').strip()
        targets  = row.get('target', '').strip()

        if not targets:
            skipped_no_target += 1
            continue

        for symbol in targets.split('|'):
            symbol = symbol.strip()
            if not symbol:
                continue
            ensg_set = symbol2ensg.get(symbol)
            if not ensg_set:
                skipped_no_mapping += 1
                continue
            for ensg in ensg_set:
                edges.add((compound, ensg))

    print(f'[repohub] Edges antes de filtro outliers: {len(edges):,}')
    print(f'[repohub] Sin target:     {skipped_no_target:,}')
    print(f'[repohub] Sin mapeo ENSG: {skipped_no_mapping:,}')

    # Eliminar outliers
    n1_out, n2_out = find_outliers_from_edges(edges)
    print(f'[repohub] Outliers n1 (compuestos): {len(n1_out):,}')
    print(f'[repohub] Outliers n2 (genes):      {len(n2_out):,}')
    edges = {(n1, n2) for n1, n2 in edges if n1 not in n1_out and n2 not in n2_out}
    print(f'[repohub] Edges tras filtro outliers: {len(edges):,}')

    # Guardar
    out_path = os.path.join(out_dir, 'CPD-int-GEN.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\n')
        for cpd, ensg in sorted(edges):
            f.write(f'{cpd}\t{ensg}\n')

    print(f'[repohub] CPD-int-GEN -> {out_path}')

    # Stats
    df_out = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':       'repohub',
        'relation':       'CPD-int-GEN',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'pert_iname (Drug Repurposing Hub)',
        'total_assoc':    len(df_out),
        'unique_genes':   df_out['n2'].nunique(),
        'unique_compounds': df_out['n1'].nunique(),
        'db_total':       len(df_out),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[repohub] Procesado completado.')


def main():
    args = parse_args()
    symbol2ensg = load_symbol2ensg(args.symbol2ensg)
    df = load_data(args.raw_dir)
    process(df, args.out_dir, symbol2ensg)
    sys.stderr.write('[repohub] Done!\n')


if __name__ == '__main__':
    main()
