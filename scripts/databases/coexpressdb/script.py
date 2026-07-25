#!/usr/bin/env python3
"""
coexpressdb -- script.py
Procesado de COXPRESdb v22-05 (Hsa-r.c6-0).

Logica (fiel a Bioteque):
  1. Cargar genes confiables desde supportability (Rank > 0, Platform Hsa-r)
  2. Para cada gen fuente, leer su archivo de coexpresion
  3. Filtrar pares cuyo mutual_rank <= cutoff (top 10% de genes totales)
  4. El target tambien debe estar en reliable_genes
  5. Mapear Entrez -> ENSG
  6. Si un par aparece varias veces (por multiples ENSG), promediar con gmean
  7. Deduplicar pares (ordenar el par para que sea simetrico)

Parametros (identicos a Bioteque):
  cutoff = ceil(n_genes * 0.1)   (~10% del total de genes)

Diferencias respecto a Bioteque:
  - Bioteque mapeaba Entrez -> UniProt; nosotros mapeamos Entrez -> ENSG
  - Supportability: Bioteque usaba supportability.2014-08-19.txt con columna Hsa3;
    nosotros usamos NonUnionPlatforms con Platform=Hsa-r y Rank > 0

Columnas archivo de coexpresion (un archivo por gen, nombre = Entrez ID):
  target_entrez_id | mutual_rank

Columnas supportability (NonUnionPlatforms):
  GeneID | Rank | Platform | Reference | MaxCoxSim

Mapeo: Entrez -> ENSG via data/processed/entrez2ensg.tsv

Output:
  processed/coexpressdb/GEN-cex-GEN.tsv   (n1=ENSG, n2=ENSG, mutual_rank)
  processed/coexpressdb/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --entrez2ensg <ruta_entrez2ensg.tsv>
"""

import argparse
import os
import sys
import math
import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',     required=True)
    parser.add_argument('--out-dir',     required=True)
    parser.add_argument('--entrez2ensg', required=True,
                        help='Ruta a data/processed/entrez2ensg.tsv')
    return parser.parse_args()


def load_entrez2ensg(path: str) -> dict:
    """Carga mapeo Entrez -> set(ENSG)."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    mapping = df.groupby('entrez_gene_id')['ensembl_gene_id'].apply(list).to_dict()
    print(f'[coexpressdb] Mapeo Entrez->ENSG: {len(mapping):,} entradas')
    return mapping


def load_reliable_genes(raw_dir: str) -> set:
    """
    Carga genes confiables desde el archivo de supportability.
    Usa NonUnionPlatforms, Platform=Hsa-r, Rank > 0.
    """
    candidates = [
        f for f in os.listdir(raw_dir)
        if 'NonUnion' in f and f.endswith('.tsv')
    ]
    if not candidates:
        raise FileNotFoundError(
            f'No se encontro archivo NonUnionPlatforms en {raw_dir}'
        )

    path = os.path.join(raw_dir, candidates[0])
    print(f'[coexpressdb] Cargando supportability desde {path} ...')

    reliable = set()
    with open(path, 'r') as f:
        next(f)  # cabecera
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            gene_id  = cols[0].strip()
            rank     = cols[1].strip()
            platform = cols[2].strip()
            if platform.startswith('Hsa-r') and rank != '0':
                reliable.add(gene_id)

    print(f'[coexpressdb] Genes confiables (Hsa-r, Rank>0): {len(reliable):,}')
    return reliable


def process(raw_dir: str, out_dir: str, entrez2ensg: dict, reliable_genes: set):
    os.makedirs(out_dir, exist_ok=True)

    coex_dir = os.path.join(raw_dir, 'coex_data')
    if not os.path.isdir(coex_dir):
        raise FileNotFoundError(
            f'No se encontro directorio coex_data en {raw_dir}. '
            f'Descomprime el zip primero.'
        )

    all_genes = os.listdir(coex_dir)
    n_genes   = len(all_genes)
    cutoff    = math.ceil(n_genes * 0.1)
    print(f'[coexpressdb] Genes totales: {n_genes:,} | Cutoff (10%): {cutoff:,}')

    # Solo procesar genes que tienen mapeo Entrez->ENSG
    genes_to_process = [g for g in all_genes if g in entrez2ensg]
    print(f'[coexpressdb] Genes con mapeo ENSG: {len(genes_to_process):,}')

    pairs = {}
    skipped_no_reliable = 0
    skipped_above_cutoff = 0

    for source_gene in tqdm(genes_to_process, desc='[coexpressdb] Procesando genes'):
        gene_file = os.path.join(coex_dir, source_gene)
        with open(gene_file, 'r') as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 2:
                    continue
                target_gene = cols[0].strip()
                try:
                    rank = float(cols[1].strip())
                except ValueError:
                    continue

                if rank < 0:
                    rank = 0.0

                if rank > cutoff:
                    skipped_above_cutoff += 1
                    continue

                if target_gene not in reliable_genes:
                    skipped_no_reliable += 1
                    continue

                if target_gene not in entrez2ensg:
                    continue

                # Mapear ambos genes a ENSG
                for ensg1 in entrez2ensg[source_gene]:
                    for ensg2 in entrez2ensg[target_gene]:
                        if ensg1 == ensg2:
                            continue
                        # Par simetrico (igual que Bioteque)
                        pair = tuple(sorted([ensg1, ensg2]))
                        if pair in pairs:
                            pairs[pair] = (pairs[pair] * rank) ** 0.5
                        else:
                            pairs[pair] = rank

    print(f'[coexpressdb] Pares unicos: {len(pairs):,}')
    print(f'[coexpressdb] Saltados (rank > cutoff): {skipped_above_cutoff:,}')
    print(f'[coexpressdb] Saltados (no reliable):   {skipped_no_reliable:,}')

    # Guardar
    out_path = os.path.join(out_dir, 'GEN-cex-GEN.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\tmutual_rank\n')
        for (ensg1, ensg2), rank in sorted(pairs.items()):
            f.write(f'{ensg1}\t{ensg2}\t{int(rank)}\n')

    print(f'[coexpressdb] GEN-cex-GEN -> {out_path}')

    # Stats
    df = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':       'coexpressdb',
        'relation':       'GEN-cex-GEN',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'ENSG',
        'total_assoc':    len(df),
        'unique_genes_n1': df['n1'].nunique(),
        'unique_genes_n2': df['n2'].nunique(),
        'cutoff_used':    cutoff,
        'db_total':       len(df),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[coexpressdb] Procesado completado.')


def main():
    args = parse_args()
    entrez2ensg    = load_entrez2ensg(args.entrez2ensg)
    reliable_genes = load_reliable_genes(args.raw_dir)
    process(args.raw_dir, args.out_dir, entrez2ensg, reliable_genes)
    sys.stderr.write('[coexpressdb] Done!\n')


if __name__ == '__main__':
    main()
