#!/usr/bin/env python3
"""
opentargets -- script.py
Procesado de Open Targets Platform 26.03 (formato Parquet).

Logica (fiel a Bioteque, adaptada a formato Parquet):
  Para cada enfermedad:
    1. Tomar todos los genes con score >= sign_cutoff (0.7)
    2. Ordenar el resto por score descendente y anadir hasta minimum_hits (100) genes
    3. Ignorar asociaciones con score < minimum_score (0.3)

Parametros (identicos a Bioteque):
  sign_cutoff   = 0.7
  minimum_hits  = 100
  minimum_score = 0.3

Estrategia de memoria (procesado por chunks):
  Se procesa un parquet a la vez. Las asociaciones de cada enfermedad se
  acumulan en un diccionario en memoria. Al final se aplica la logica
  de seleccion y se escribe el output.
  Pico de memoria estimado: proporcional al numero de asociaciones unicas
  por enfermedad, no al tamano total del dataset.

Columnas de los parquets:
  diseaseId | targetId | associationScore

Normalizacion de IDs:
  diseaseId: DOID_0050890  -> DOID:0050890
             Orphanet_166  -> ORPHA:166
  targetId:  ENSG (sin mapeo necesario)

Output:
  processed/opentargets/GEN-ass-DIS.tsv   (n1=ENSG, n2=DiseaseID, score)
  processed/opentargets/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
"""

import argparse
import os
import sys
import glob
import pandas as pd
from collections import defaultdict


# Parametros identicos a Bioteque
SIGN_CUTOFF   = 0.7
MINIMUM_HITS  = 100
MINIMUM_SCORE = 0.3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    return parser.parse_args()


def normalize_disease_id(dis_id: str) -> str:
    """
    Normaliza el ID de enfermedad de Open Targets.
      DOID_0050890 -> DOID:0050890
      EFO_0000270  -> EFO:0000270
      Orphanet_166 -> ORPHA:166
      HP_0001250   -> HP:0001250
    """
    if dis_id.startswith('Orphanet_'):
        return 'ORPHA:' + dis_id.split('_', 1)[1]
    return dis_id.replace('_', ':', 1)


def load_associations_chunked(raw_dir: str) -> dict:
    """
    Lee los parquets de uno en uno y acumula asociaciones por enfermedad.
    Retorna: {disease_id: {(ensg, score), ...}}
    Solo guarda pares que superen MINIMUM_SCORE para ahorrar memoria.
    """
    parquet_files = sorted(glob.glob(os.path.join(raw_dir, '*.parquet')))
    if not parquet_files:
        raise FileNotFoundError(f'No se encontraron parquets en {raw_dir}')

    print(f'[opentargets] Procesando {len(parquet_files)} parquets en chunks...')

    # dis_id -> dict de ensg -> score (guardamos el score mas alto si hay duplicados)
    dis_ass = defaultdict(dict)

    for i, fpath in enumerate(parquet_files):
        print(f'  [{i+1}/{len(parquet_files)}] {os.path.basename(fpath)}')
        df = pd.read_parquet(fpath, columns=['diseaseId', 'targetId', 'associationScore'])

        # Filtro rapido antes de acumular
        df = df[df['associationScore'] >= MINIMUM_SCORE]

        for _, row in df.iterrows():
            dis = normalize_disease_id(row['diseaseId'])
            ensg  = row['targetId']
            score = row['associationScore']

            # Conservar score mas alto si el par ya existe
            if ensg not in dis_ass[dis] or dis_ass[dis][ensg] < score:
                dis_ass[dis][ensg] = score

    print(f'[opentargets] Enfermedades acumuladas: {len(dis_ass):,}')
    return dis_ass


def write_output(dis_ass: dict, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'GEN-ass-DIS.tsv')

    print('[opentargets] Aplicando logica sign_cutoff / minimum_hits y escribiendo...')

    total_written = 0
    unique_genes  = set()
    unique_dis    = set()

    with open(out_path, 'w') as o:
        o.write('n1\tn2\tscore\n')

        for dis_id, gene_scores in dis_ass.items():
            # Ordenar genes por score descendente
            genes_sorted = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)

            c = 0
            for ensg, score in genes_sorted:
                if score < MINIMUM_SCORE:
                    break
                elif score >= SIGN_CUTOFF:
                    o.write(f'{ensg}\t{dis_id}\t{score:.3f}\n')
                    unique_genes.add(ensg)
                    unique_dis.add(dis_id)
                    total_written += 1
                    c += 1
                else:
                    if c < MINIMUM_HITS:
                        o.write(f'{ensg}\t{dis_id}\t{score:.3f}\n')
                        unique_genes.add(ensg)
                        unique_dis.add(dis_id)
                        total_written += 1
                        c += 1
                    else:
                        break

    print(f'[opentargets] GEN-ass-DIS: {total_written:,} asociaciones -> {out_path}')
    print(f'[opentargets] Genes unicos:        {len(unique_genes):,}')
    print(f'[opentargets] Enfermedades unicas: {len(unique_dis):,}')

    # Stats
    stats = pd.DataFrame([{
        'database':       'opentargets',
        'relation':       'GEN-ass-DIS',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'EFO/DOID/ORPHA/HP',
        'total_assoc':    total_written,
        'unique_genes':   len(unique_genes),
        'unique_targets': len(unique_dis),
        'db_total':       total_written,
        'sign_cutoff':    SIGN_CUTOFF,
        'minimum_hits':   MINIMUM_HITS,
        'minimum_score':  MINIMUM_SCORE,
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[opentargets] Procesado completado.')


def main():
    args = parse_args()
    dis_ass = load_associations_chunked(args.raw_dir)
    write_output(dis_ass, args.out_dir)
    sys.stderr.write('[opentargets] Done!\n')


if __name__ == '__main__':
    main()
