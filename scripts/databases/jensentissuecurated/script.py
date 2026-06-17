#!/usr/bin/env python3
"""
jensentissuecurated -- script.py
Procesado de TISSUES Jensen Lab (human_tissue_knowledge_full.tsv).

Logica (fiel a Bioteque):
  1. Leer archivo de asociaciones gen-tejido
  2. Saltar BTO:0000000 (raiz de la ontologia)
  3. Mapear ENSP -> ENSG via archivo de Ensembl
  4. Deduplicar pares (ENSG, BTO_ID)

Columnas del archivo fuente:
  ENSP | gene_symbol | BTO_ID | tissue_name | source | evidence | score

Mapeo: ENSP -> ENSG via data/processed/ensp2ensg.tsv.gz
  (Homo_sapiens.GRCh38.113.uniprot.tsv.gz de Ensembl FTP)
  Columnas: gene_stable_id | transcript_stable_id | protein_stable_id | ...

Output:
  processed/jensentissuecurated/GEN-ass-TIS.tsv   (n1=ENSG, n2=BTO_ID)
  processed/jensentissuecurated/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --ensp2ensg <ruta_ensp2ensg.tsv.gz>
"""

import argparse
import gzip
import os
import sys
import pandas as pd


BTO_ROOT = 'BTO:0000000'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',   required=True)
    parser.add_argument('--out-dir',   required=True)
    parser.add_argument('--ensp2ensg', required=True,
                        help='Ruta a data/processed/ensp2ensg.tsv.gz')
    return parser.parse_args()


def load_ensp2ensg(path: str) -> dict:
    """
    Carga mapeo ENSP -> ENSG desde el archivo de Ensembl.
    Un ENSP puede mapear a un solo ENSG.
    """
    print(f'[jensentissuecurated] Cargando mapeo ENSP->ENSG desde {path} ...')
    mapping = {}
    with gzip.open(path, 'rt') as f:
        next(f)  # saltar cabecera
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            ensg = cols[0].strip()
            ensp = cols[2].strip()
            if ensg.startswith('ENSG') and ensp.startswith('ENSP'):
                if ensp not in mapping:
                    mapping[ensp] = ensg
    print(f'[jensentissuecurated] Mapeo cargado: {len(mapping):,} ENSP -> ENSG')
    return mapping


def process(raw_dir: str, out_dir: str, ensp2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)
    infile = os.path.join(raw_dir, 'human_tissue_knowledge_full.tsv')

    print(f'[jensentissuecurated] Procesando {infile} ...')

    pairs = set()
    skipped_root    = 0
    skipped_mapping = 0
    n_lines         = 0

    with open(infile, 'r') as f:
        for line in f:
            n_lines += 1
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue

            ensp   = cols[0].strip()
            bto_id = cols[2].strip()

            # Saltar raiz de ontologia (igual que Bioteque)
            if bto_id == BTO_ROOT:
                skipped_root += 1
                continue

            # Mapear ENSP -> ENSG
            ensg = ensp2ensg.get(ensp)
            if ensg is None:
                skipped_mapping += 1
                continue

            pairs.add((ensg, bto_id))

    print(f'[jensentissuecurated] Lineas procesadas:      {n_lines:,}')
    print(f'[jensentissuecurated] Saltadas (BTO root):    {skipped_root:,}')
    print(f'[jensentissuecurated] Sin mapeo ENSP->ENSG:   {skipped_mapping:,}')
    print(f'[jensentissuecurated] Pares unicos (ENSG,BTO): {len(pairs):,}')

    # Guardar
    out_path = os.path.join(out_dir, 'GEN-ass-TIS.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, bto in sorted(pairs):
            f.write(f'{ensg}\t{bto}\n')

    print(f'[jensentissuecurated] GEN-ass-TIS -> {out_path}')

    # Stats
    df = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':       'jensentissuecurated',
        'relation':       'GEN-ass-TIS',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'BTO ID',
        'total_assoc':    len(df),
        'unique_genes':   df['n1'].nunique(),
        'unique_targets': df['n2'].nunique(),
        'db_total':       len(df),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[jensentissuecurated] Procesado completado.')


def main():
    args = parse_args()
    ensp2ensg = load_ensp2ensg(args.ensp2ensg)
    process(args.raw_dir, args.out_dir, ensp2ensg)
    sys.stderr.write('[jensentissuecurated] Done!\n')


if __name__ == '__main__':
    main()
