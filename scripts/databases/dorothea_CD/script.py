#!/usr/bin/env python3
"""
dorothea_CD -- script.py
Procesado de DoRothEA niveles C y D via OmniPath interactions.

Logica (fiel a Bioteque):
  1. Filtrar df[dorothea == True]
  2. Filtrar solo humano (ncbi_tax_id_source == 9606 y ncbi_tax_id_target == 9606)
  3. Filtrar source y target en UniProt revisados
  4. Filtrar niveles de confianza A y B
  5. Separar en reg/upr/dwr por consensus_stimulation y consensus_inhibition

Mapeo: UniProt -> ENSG via ensp2ensg.tsv.gz (Uniprot/SWISSPROT)

Niveles de confianza DoRothEA:
  A: ChIP-seq experimental + literatura curada
  B: ChIP-seq o literatura curada

Output:
  processed/dorothea_CD/GEN-reg-GEN.tsv
  processed/dorothea_CD/GEN-upr-GEN.tsv
  processed/dorothea_CD/GEN-dwr-GEN.tsv
  processed/dorothea_CD/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --ensp2ensg <ruta_ensp2ensg.tsv.gz>
"""

import argparse
import gzip
import os
import sys
import pandas as pd


DOROTHEA_LEVELS = ["C", "D"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',   required=True)
    parser.add_argument('--out-dir',   required=True)
    parser.add_argument('--ensp2ensg', required=True)
    return parser.parse_args()


def load_reviewed_uniprot(path: str) -> set:
    reviewed = set()
    with gzip.open(path, 'rt') as f:
        next(f)
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 5 and cols[4].strip() == 'Uniprot/SWISSPROT':
                reviewed.add(cols[3].strip())
    print(f'[dorothea_CD] UniProt revisados: {len(reviewed):,}')
    return reviewed


def load_uniprot2ensg(path: str) -> dict:
    mapping = {}
    with gzip.open(path, 'rt') as f:
        next(f)
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue
            ensg    = cols[0].strip()
            uniprot = cols[3].strip()
            db_name = cols[4].strip()
            if db_name != 'Uniprot/SWISSPROT' or not ensg.startswith('ENSG'):
                continue
            if uniprot not in mapping:
                mapping[uniprot] = set()
            mapping[uniprot].add(ensg)
    print(f'[dorothea_CD] Mapeo UniProt->ENSG: {len(mapping):,} entradas')
    return mapping


def to_ensg_pairs(df: pd.DataFrame, up2ensg: dict) -> set:
    pairs = set()
    for _, row in df.iterrows():
        for e1 in up2ensg.get(row['source'], set()):
            for e2 in up2ensg.get(row['target'], set()):
                if e1 != e2:
                    pairs.add((e1, e2))
    return pairs


def process(raw_dir: str, out_dir: str, reviewed: set, up2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)

    infile = os.path.join(raw_dir, 'omnipath_webservice_interactions__latest.tsv.gz')
    print(f'[dorothea_CD] Cargando {infile} ...')

    df = pd.read_csv(infile, sep='\t', low_memory=False)
    print(f'  Filas totales: {len(df):,}')

    # Filtrar dorothea == True
    df = df[df['dorothea'] == True]
    print(f'  Tras filtro dorothea=True: {len(df):,}')

    # Filtrar humano
    df = df[(df['ncbi_tax_id_source'] == 9606) & (df['ncbi_tax_id_target'] == 9606)]
    print(f'  Tras filtro humano: {len(df):,}')

    # Filtrar revisados
    df = df[(df['source'].isin(reviewed)) & (df['target'].isin(reviewed))]
    print(f'  Tras filtro revisados: {len(df):,}')

    # Filtrar niveles C y D (igual que Bioteque: busca si el nivel esta en el string)
    levels = set()
    for x in df['dorothea_level'].dropna().unique():
        for lvl in DOROTHEA_LEVELS:
            if lvl in x:
                levels.add(x)
    df = df[df['dorothea_level'].isin(levels)]
    print(f'  Tras filtro niveles {DOROTHEA_LEVELS}: {len(df):,}')

    # Regulacion general
    reg = to_ensg_pairs(df, up2ensg)

    # Up-regulacion
    upr = to_ensg_pairs(df[df['consensus_stimulation'] == True], up2ensg)

    # Down-regulacion
    dwr = to_ensg_pairs(df[df['consensus_inhibition'] == True], up2ensg)

    # Guardar
    for name, pairs in [('GEN-reg-GEN', reg), ('GEN-upr-GEN', upr), ('GEN-dwr-GEN', dwr)]:
        path = os.path.join(out_dir, f'{name}.tsv')
        with open(path, 'w') as f:
            f.write('n1\tn2\n')
            for e1, e2 in sorted(pairs):
                f.write(f'{e1}\t{e2}\n')
        print(f'[dorothea_CD] {name}: {len(pairs):,} -> {path}')

    # Stats
    total_db = len(reg) + len(upr) + len(dwr)
    stats_rows = []
    for rel, s in [('GEN-reg-GEN', reg), ('GEN-upr-GEN', upr), ('GEN-dwr-GEN', dwr)]:
        stats_rows.append({
            'database':      'dorothea_CD',
            'relation':      rel,
            'gene_id_type':  'ENSG',
            'assoc_id_type': 'ENSG',
            'total_assoc':   len(s),
            'unique_genes':  len({e for pair in s for e in pair}),
            'db_total':      total_db,
        })
    pd.DataFrame(stats_rows).to_csv(
        os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False
    )
    print('[dorothea_CD] Procesado completado.')


def main():
    args = parse_args()
    reviewed  = load_reviewed_uniprot(args.ensp2ensg)
    up2ensg   = load_uniprot2ensg(args.ensp2ensg)
    process(args.raw_dir, args.out_dir, reviewed, up2ensg)
    sys.stderr.write('[dorothea_CD] Done!\n')


if __name__ == '__main__':
    main()
