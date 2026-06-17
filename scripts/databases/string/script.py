#!/usr/bin/env python3
"""
string -- script.py
Procesado de STRING v12.0 (9606.protein.links.full.v12.0.txt.gz).

Logica (fiel a Bioteque):
  1. Filtrar pares con combined_score >= 700 (equivalente a 0.7)
  2. Mapear ENSP -> ENSG via ensp2ensg.tsv.gz
  3. Pares simetricos (ordenar el par)
  4. Si un par aparece varias veces, guardar el score minimo (igual que Bioteque)

Diferencias respecto a Bioteque:
  - Bioteque mapeaba STRING ID -> UniProt -> UniProt revisados
  - Nosotros mapeamos ENSP -> ENSG directamente (mas simple y actualizado)
  - El archivo nuevo usa ENSP en vez de STRING protein IDs

Parametros:
  min_score = 0.7  (combined_score >= 700 en escala 0-1000)

Columnas del archivo:
  protein1 | protein2 | neighborhood | ... | combined_score
  (protein1/2 formato: 9606.ENSP00000000233)

Output:
  processed/string/GEN-ppi-GEN.tsv   (n1=ENSG, n2=ENSG, score)
  processed/string/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --ensp2ensg <ruta_ensp2ensg.tsv.gz>
"""

import argparse
import gzip
import os
import sys
import pandas as pd
from tqdm import tqdm


MIN_SCORE = 700  # equivalente a 0.7 en escala 0-1000


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',   required=True)
    parser.add_argument('--out-dir',   required=True)
    parser.add_argument('--ensp2ensg', required=True,
                        help='Ruta a data/processed/ensp2ensg.tsv.gz')
    return parser.parse_args()


def load_ensp2ensg(path: str) -> dict:
    """Carga mapeo ENSP -> ENSG (solo un ENSG por ENSP)."""
    print(f'[string] Cargando mapeo ENSP->ENSG desde {path} ...')
    mapping = {}
    with gzip.open(path, 'rt') as f:
        next(f)
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            ensg = cols[0].strip()
            ensp = cols[2].strip()
            if ensg.startswith('ENSG') and ensp.startswith('ENSP'):
                if ensp not in mapping:
                    mapping[ensp] = ensg
    print(f'[string] Mapeo ENSP->ENSG: {len(mapping):,} entradas')
    return mapping


def process(raw_dir: str, out_dir: str, ensp2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)

    infile = os.path.join(raw_dir, '9606.protein.links.full.v12.0.txt.gz')
    print(f'[string] Procesando {infile} ...')

    edges = {}
    skipped_score   = 0
    skipped_mapping = 0
    n_lines         = 0

    with gzip.open(infile, 'rt') as f:
        f.readline()  # cabecera

        for line in tqdm(f, desc='[string] Leyendo pares'):
            n_lines += 1
            cols = line.rstrip('\n').split(' ')
            if len(cols) < 16:
                continue

            score = int(cols[-1])
            if score < MIN_SCORE:
                skipped_score += 1
                continue

            # Extraer ENSP (quitar prefijo 9606.)
            ensp1 = cols[0].split('.')[-1]
            ensp2 = cols[1].split('.')[-1]

            ensg1 = ensp2ensg.get(ensp1)
            ensg2 = ensp2ensg.get(ensp2)

            if ensg1 is None or ensg2 is None:
                skipped_mapping += 1
                continue

            if ensg1 == ensg2:
                continue

            pair = tuple(sorted([ensg1, ensg2]))
            if pair not in edges:
                edges[pair] = []
            edges[pair].append(score / 1000)

    print(f'[string] Lineas procesadas:        {n_lines:,}')
    print(f'[string] Saltadas (score < 0.7):   {skipped_score:,}')
    print(f'[string] Saltadas (sin mapeo ENSG): {skipped_mapping:,}')
    print(f'[string] Pares unicos:             {len(edges):,}')

    # Guardar con score minimo (igual que Bioteque)
    out_path = os.path.join(out_dir, 'GEN-ppi-GEN.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\tscore\n')
        for (ensg1, ensg2), scores in sorted(edges.items()):
            f.write(f'{ensg1}\t{ensg2}\t{min(scores):.3f}\n')

    print(f'[string] GEN-ppi-GEN -> {out_path}')

    # Stats
    df = pd.read_csv(out_path, sep='\t')
    all_genes = set(df['n1']) | set(df['n2'])
    stats = pd.DataFrame([{
        'database':      'string',
        'relation':      'GEN-ppi-GEN',
        'gene_id_type':  'ENSG',
        'assoc_id_type': 'ENSG',
        'total_assoc':   len(df),
        'unique_genes':  len(all_genes),
        'min_score':     MIN_SCORE / 1000,
        'db_total':      len(df),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[string] Procesado completado.')


def main():
    args = parse_args()
    ensp2ensg = load_ensp2ensg(args.ensp2ensg)
    process(args.raw_dir, args.out_dir, ensp2ensg)
    sys.stderr.write('[string] Done!\n')


if __name__ == '__main__':
    main()
