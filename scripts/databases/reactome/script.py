#!/usr/bin/env python3
"""
reactome -- script.py
Procesado de Reactome (UniProt2Reactome_All_Levels.txt).

Logica (fiel a Bioteque):
  1. Filtrar solo pathways humanos (R-HSA-)
  2. Filtrar solo UniProt IDs revisados (Uniprot/SWISSPROT)
  3. Mapear UniProt -> ENSG via archivo de Ensembl
  4. Deduplicar pares (ENSG, ReactomeID)

Columnas del archivo fuente:
  UniProtID | ReactomeID | URL | pathway_name | evidence | species

Mapeo: UniProt -> ENSG via data/processed/ensp2ensg.tsv.gz
  Solo se usan entradas con db_name == Uniprot/SWISSPROT (revisados)
  para ser consistente con la logica de Bioteque (human reviewed UniProt).

Output:
  processed/reactome/GEN-ass-PWY.tsv   (n1=ENSG, n2=ReactomeID)
  processed/reactome/stats.tsv

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


def load_uniprot2ensg(path: str) -> dict:
    """
    Carga mapeo UniProt -> ENSG desde el archivo de Ensembl.
    Solo usa entradas Uniprot/SWISSPROT (revisadas, igual que Bioteque).
    Un UniProt puede mapear a varios ENSG.
    """
    print(f'[reactome] Cargando mapeo UniProt->ENSG desde {path} ...')
    mapping = {}
    with gzip.open(path, 'rt') as f:
        next(f)  # saltar cabecera
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue
            ensg    = cols[0].strip()
            uniprot = cols[3].strip()
            db_name = cols[4].strip()

            # Solo UniProt revisados (SWISSPROT), igual que Bioteque
            if db_name != 'Uniprot/SWISSPROT':
                continue
            if not ensg.startswith('ENSG') or not uniprot:
                continue

            if uniprot not in mapping:
                mapping[uniprot] = set()
            mapping[uniprot].add(ensg)

    print(f'[reactome] Mapeo cargado: {len(mapping):,} UniProt -> ENSG')
    return mapping


def process(raw_dir: str, out_dir: str, uniprot2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)
    infile = os.path.join(raw_dir, 'UniProt2Reactome_All_Levels.txt')

    print(f'[reactome] Procesando {infile} ...')

    pairs = set()
    skipped_species = 0
    skipped_mapping = 0
    n_lines         = 0

    with open(infile, 'r') as f:
        for line in f:
            n_lines += 1
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 2:
                continue

            uniprot    = cols[0].strip()
            reactome_id = cols[1].strip()

            # Solo pathways humanos (igual que Bioteque: startswith R-HSA-)
            if not reactome_id.startswith('R-HSA-'):
                skipped_species += 1
                continue

            # Mapear UniProt -> ENSG
            ensg_set = uniprot2ensg.get(uniprot)
            if not ensg_set:
                skipped_mapping += 1
                continue

            for ensg in ensg_set:
                pairs.add((ensg, reactome_id))

    print(f'[reactome] Lineas procesadas:        {n_lines:,}')
    print(f'[reactome] Saltadas (no humano):     {skipped_species:,}')
    print(f'[reactome] Sin mapeo UniProt->ENSG:  {skipped_mapping:,}')
    print(f'[reactome] Pares unicos (ENSG, PWY): {len(pairs):,}')

    # Guardar
    out_path = os.path.join(out_dir, 'GEN-ass-PWY.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, pwy in sorted(pairs, key=lambda x: x[1]):
            f.write(f'{ensg}\t{pwy}\n')

    print(f'[reactome] GEN-ass-PWY -> {out_path}')

    # Stats
    df = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':       'reactome',
        'relation':       'GEN-ass-PWY',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'Reactome Pathway ID (R-HSA-)',
        'total_assoc':    len(df),
        'unique_genes':   df['n1'].nunique(),
        'unique_targets': df['n2'].nunique(),
        'db_total':       len(df),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[reactome] Procesado completado.')


def main():
    args = parse_args()
    uniprot2ensg = load_uniprot2ensg(args.ensp2ensg)
    process(args.raw_dir, args.out_dir, uniprot2ensg)
    sys.stderr.write('[reactome] Done!\n')


if __name__ == '__main__':
    main()
