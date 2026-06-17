#!/usr/bin/env python3
"""
omnipath -- script.py
Procesado de OmniPath interactions y enzyme-substrate.

Logica (fiel a Bioteque):
  PPIs:
    1. Filtrar df[omnipath == True]
    2. Filtrar solo humano (ncbi_tax_id_source == 9606 y ncbi_tax_id_target == 9606)
    3. Filtrar source y target en UniProt revisados (Uniprot/SWISSPROT)
    4. Pares simetricos (ordenar el par)

  PTMs (enzyme-substrate):
    1. Filtrar solo humano (ncbi_tax_id == 9606)
    2. Filtrar enzyme y substrate en UniProt revisados
    3. Separar por modification: phosphorylation -> GEN-pho-GEN
                                  dephosphorylation -> GEN-dph-GEN

  Mapeo UniProt -> ENSG via ensp2ensg.tsv.gz (solo Uniprot/SWISSPROT)

Columnas interactions:
  source | target | source_genesymbol | target_genesymbol | is_directed |
  ... | omnipath | ncbi_tax_id_source | ncbi_tax_id_target | ...

Columnas enz_sub:
  enzyme | substrate | ... | ncbi_tax_id | modification | ...

Output:
  processed/omnipath/GEN-ppi-GEN.tsv
  processed/omnipath/GEN-pho-GEN.tsv
  processed/omnipath/GEN-dph-GEN.tsv
  processed/omnipath/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --ensp2ensg <ruta_ensp2ensg.tsv.gz>
"""

import argparse
import gzip
import os
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',   required=True)
    parser.add_argument('--out-dir',   required=True)
    parser.add_argument('--ensp2ensg', required=True,
                        help='Ruta a data/processed/ensp2ensg.tsv.gz')
    return parser.parse_args()


def load_reviewed_uniprot(path: str) -> set:
    """Carga UniProt IDs revisados (Uniprot/SWISSPROT) desde ensp2ensg.tsv.gz."""
    print(f'[omnipath] Cargando UniProt revisados desde {path} ...')
    reviewed = set()
    with gzip.open(path, 'rt') as f:
        next(f)
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue
            if cols[4].strip() == 'Uniprot/SWISSPROT':
                reviewed.add(cols[3].strip())
    print(f'[omnipath] UniProt revisados: {len(reviewed):,}')
    return reviewed


def load_uniprot2ensg(path: str) -> dict:
    """Carga mapeo UniProt -> set(ENSG) solo para revisados."""
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
            if db_name != 'Uniprot/SWISSPROT':
                continue
            if not ensg.startswith('ENSG'):
                continue
            if uniprot not in mapping:
                mapping[uniprot] = set()
            mapping[uniprot].add(ensg)
    print(f'[omnipath] Mapeo UniProt->ENSG: {len(mapping):,} entradas')
    return mapping


def write_pairs(pairs: set, path: str, col1: str = 'n1', col2: str = 'n2'):
    with open(path, 'w') as f:
        f.write(f'{col1}\t{col2}\n')
        for a, b in sorted(pairs):
            f.write(f'{a}\t{b}\n')


def process_ppi(raw_dir: str, reviewed: set, up2ensg: dict) -> set:
    """Procesa interactions -> GEN-ppi-GEN."""
    infile = os.path.join(raw_dir, 'omnipath_webservice_interactions__latest.tsv.gz')
    print(f'[omnipath] Cargando interactions desde {infile} ...')

    df = pd.read_csv(infile, sep='\t', low_memory=False)
    print(f'  Filas totales: {len(df):,}')

    # Verificar columnas de tax_id
    tax_cols = [c for c in df.columns if 'tax' in c.lower()]
    print(f'  Columnas tax: {tax_cols}')

    # Filtrar omnipath == True
    df = df[df['omnipath'] == True]
    print(f'  Tras filtro omnipath=True: {len(df):,}')

    # Filtrar humano
    if 'ncbi_tax_id_source' in df.columns and 'ncbi_tax_id_target' in df.columns:
        df = df[(df['ncbi_tax_id_source'] == 9606) & (df['ncbi_tax_id_target'] == 9606)]
    print(f'  Tras filtro humano: {len(df):,}')

    # Filtrar UniProt revisados
    df = df[(df['source'].isin(reviewed)) & (df['target'].isin(reviewed))]
    print(f'  Tras filtro revisados: {len(df):,}')

    # Mapear a ENSG y construir pares simetricos
    ppis = set()
    for _, row in df.iterrows():
        for e1 in up2ensg.get(row['source'], set()):
            for e2 in up2ensg.get(row['target'], set()):
                if e1 != e2:
                    ppis.add(tuple(sorted([e1, e2])))

    print(f'[omnipath] GEN-ppi-GEN: {len(ppis):,} pares')
    return ppis


def process_ptm(raw_dir: str, reviewed: set, up2ensg: dict) -> tuple:
    """Procesa enzyme-substrate -> GEN-pho-GEN y GEN-dph-GEN."""
    infile = os.path.join(raw_dir, 'omnipath_webservice_enz_sub__latest.tsv.gz')
    print(f'[omnipath] Cargando enz_sub desde {infile} ...')

    df = pd.read_csv(infile, sep='\t', low_memory=False)
    print(f'  Filas totales: {len(df):,}')

    # Filtrar humano
    if 'ncbi_tax_id' in df.columns:
        df = df[df['ncbi_tax_id'] == 9606]
    print(f'  Tras filtro humano: {len(df):,}')

    # Filtrar revisados
    df = df[(df['enzyme'].isin(reviewed)) & (df['substrate'].isin(reviewed))]
    print(f'  Tras filtro revisados: {len(df):,}')

    # Separar por modificacion
    pho_df = df[df['modification'] == 'phosphorylation']
    dph_df = df[df['modification'] == 'dephosphorylation']

    def to_ensg_pairs(sub_df, symmetric=False):
        pairs = set()
        for _, row in sub_df.iterrows():
            for e1 in up2ensg.get(row['enzyme'], set()):
                for e2 in up2ensg.get(row['substrate'], set()):
                    if e1 != e2:
                        if symmetric:
                            pairs.add(tuple(sorted([e1, e2])))
                        else:
                            pairs.add((e1, e2))
        return pairs

    pho = to_ensg_pairs(pho_df)
    dph = to_ensg_pairs(dph_df)

    print(f'[omnipath] GEN-pho-GEN: {len(pho):,} pares')
    print(f'[omnipath] GEN-dph-GEN: {len(dph):,} pares')
    return pho, dph


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    reviewed  = load_reviewed_uniprot(args.ensp2ensg)
    up2ensg   = load_uniprot2ensg(args.ensp2ensg)

    # PPIs
    ppis = process_ppi(args.raw_dir, reviewed, up2ensg)
    write_pairs(ppis, os.path.join(args.out_dir, 'GEN-ppi-GEN.tsv'))

    # PTMs
    pho, dph = process_ptm(args.raw_dir, reviewed, up2ensg)
    write_pairs(pho, os.path.join(args.out_dir, 'GEN-pho-GEN.tsv'))
    write_pairs(dph, os.path.join(args.out_dir, 'GEN-dph-GEN.tsv'))

    # Stats
    total_db = len(ppis) + len(pho) + len(dph)
    stats_rows = []
    for rel, s in [('GEN-ppi-GEN', ppis), ('GEN-pho-GEN', pho), ('GEN-dph-GEN', dph)]:
        ensg_all = {e for pair in s for e in pair}
        stats_rows.append({
            'database':      'omnipath',
            'relation':      rel,
            'gene_id_type':  'ENSG',
            'assoc_id_type': 'ENSG',
            'total_assoc':   len(s),
            'unique_genes':  len(ensg_all),
            'db_total':      total_db,
        })
    pd.DataFrame(stats_rows).to_csv(
        os.path.join(args.out_dir, 'stats.tsv'), sep='\t', index=False
    )
    print('[omnipath] Procesado completado.')
    sys.stderr.write('[omnipath] Done!\n')


if __name__ == '__main__':
    main()
