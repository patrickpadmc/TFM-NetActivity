#!/usr/bin/env python3
"""
coexpressdb -- script_fix.py
Procesado de COXPRESdb v22-05 (Hsa-r.c6-0).

Fix: manejo de NaN en gmean cuando mutual_rank tiene valores invalidos.
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
    parser.add_argument('--entrez2ensg', required=True)
    return parser.parse_args()


def load_entrez2ensg(path: str) -> dict:
    df = pd.read_csv(path, sep='\t', dtype=str)
    mapping = df.groupby('entrez_gene_id')['ensembl_gene_id'].apply(list).to_dict()
    print(f'[coexpressdb] Mapeo Entrez->ENSG: {len(mapping):,} entradas')
    return mapping


def load_reliable_genes(raw_dir: str) -> set:
    candidates = [f for f in os.listdir(raw_dir) if 'NonUnion' in f and f.endswith('.tsv')]
    if not candidates:
        raise FileNotFoundError(f'No se encontro archivo NonUnionPlatforms en {raw_dir}')
    path = os.path.join(raw_dir, candidates[0])
    print(f'[coexpressdb] Cargando supportability desde {path} ...')
    reliable = set()
    with open(path, 'r') as f:
        next(f)
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
        raise FileNotFoundError(f'No se encontro directorio coex_data en {raw_dir}.')

    all_genes = os.listdir(coex_dir)
    n_genes   = len(all_genes)
    cutoff    = math.ceil(n_genes * 0.1)
    print(f'[coexpressdb] Genes totales: {n_genes:,} | Cutoff (10%): {cutoff:,}')

    genes_to_process = [g for g in all_genes if g in entrez2ensg]
    print(f'[coexpressdb] Genes con mapeo ENSG: {len(genes_to_process):,}')

    pairs = {}
    skipped_no_reliable = 0
    skipped_above_cutoff = 0
    skipped_nan = 0

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

                # Fix: saltar NaN explicitamente
                if np.isnan(rank):
                    skipped_nan += 1
                    continue

                if rank > cutoff:
                    skipped_above_cutoff += 1
                    continue

                if target_gene not in reliable_genes:
                    skipped_no_reliable += 1
                    continue

                if target_gene not in entrez2ensg:
                    continue

                for ensg1 in entrez2ensg[source_gene]:
                    for ensg2 in entrez2ensg[target_gene]:
                        if ensg1 == ensg2:
                            continue
                        pair = tuple(sorted([ensg1, ensg2]))
                        if pair in pairs:
                            new_rank = float(gmean([pairs[pair], rank]))
                            # Fix: saltar si gmean produce NaN
                            if not np.isnan(new_rank):
                                pairs[pair] = new_rank
                        else:
                            pairs[pair] = rank

    print(f'[coexpressdb] Pares unicos: {len(pairs):,}')
    print(f'[coexpressdb] Saltados (rank > cutoff): {skipped_above_cutoff:,}')
    print(f'[coexpressdb] Saltados (no reliable):   {skipped_no_reliable:,}')
    print(f'[coexpressdb] Saltados (NaN rank):       {skipped_nan:,}')

    out_path = os.path.join(out_dir, 'GEN-cex-GEN.tsv')
    written = 0
    with open(out_path, 'w') as f:
        f.write('n1\tn2\tmutual_rank\n')
        for (ensg1, ensg2), rank in sorted(pairs.items()):
            if np.isnan(rank):
                continue
            f.write(f'{ensg1}\t{ensg2}\t{int(round(rank))}\n')
            written += 1

    print(f'[coexpressdb] GEN-cex-GEN: {written:,} pares -> {out_path}')

    df = pd.read_csv(out_path, sep='\t')
    all_genes_out = set(df['n1']) | set(df['n2'])
    stats = pd.DataFrame([{
        'database':        'coexpressdb',
        'relation':        'GEN-cex-GEN',
        'gene_id_type':    'ENSG',
        'assoc_id_type':   'ENSG',
        'total_assoc':     len(df),
        'unique_genes_n1': df['n1'].nunique(),
        'unique_genes_n2': df['n2'].nunique(),
        'unique_genes':    len(all_genes_out),
        'cutoff_used':     cutoff,
        'db_total':        len(df),
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
