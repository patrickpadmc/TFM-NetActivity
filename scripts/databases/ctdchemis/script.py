#!/usr/bin/env python3
"""
ctdchemis - script.py
Procesado de CTD chemicals-diseases (CTD_chemicals_diseases.tsv.gz).

Logica (fiel a Bioteque):
  1. Mapear ChemicalID (CTD) -> InChIKey
     Bioteque usa ctd.tsv de su mapping_folder. Nosotros usamos el mismo
     archivo si esta disponible, o mantenemos el ID CTD como fallback.
  2. Filtrar por DirectEvidence:
       "marker/mechanism" -> CPD-cau-DIS
       "therapeutic"      -> CPD-trt-DIS
       (vacio/inferido) -> descartar
  3. Normalizar DiseaseID: estandarizar prefijo (MESH:, OMIM:, etc.)
     CTD usa casi exclusivamente MeSH, asi que el resultado es MESH:XXXXXXX

Columnas del archivo (linea 28 del gzip):
  ChemicalName | ChemicalID | CasRN | DiseaseName | DiseaseID |
  DirectEvidence | InferenceGeneSymbol | InferenceScore | OmimIDs | PubMedIDs

Output:
  processed/ctdchemis/CPD-cau-DIS.tsv   (n1=ChemicalID o InChIKey, n2=DiseaseID)
  processed/ctdchemis/CPD-trt-DIS.tsv
  processed/ctdchemis/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
             [--ctd-mapping <ruta_a_ctd.tsv>]
"""

import argparse
import gzip
import os
import sys
import pandas as pd


CAUSES_EV = {'marker/mechanism', 'marker', 'mechanism'}
TREATS_EV = {'therapeutic'}

# Ontologias a ignorar (igual que Bioteque)
IGNORE_SOURCES = {'MEDGEN', 'SNOMED', 'BFO', 'ICD', 'NCI'}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',     required=True)
    parser.add_argument('--out-dir',     required=True)
    parser.add_argument('--ctd-mapping', default=None,
                        help='Ruta al archivo ctd.tsv de Bioteque (col0=CTD_ID, col2=InChIKey). '
                             'Si no se pasa, se usa el ChemicalID de CTD directamente.')
    return parser.parse_args()


def load_ctd_mapping(path: str) -> dict:
    """Carga mapeo CTD ID -> InChIKey desde el archivo de Bioteque."""
    d = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3 or not cols[2]:
                continue
            d[cols[0]] = cols[2]
    print(f'[ctdchemis] Mapeo CTD -> InChIKey cargado: {len(d):,} entradas')
    return d


def normalize_disease_id(dis: str) -> str | None:
    """
    Replica la logica de parse_diseaseID de Bioteque para normalizar
    el prefijo del ID de enfermedad.
    CTD usa MESH:XXXXXXX directamente, asi que en la practica esto
    solo limpia el string.
    """
    if not dis or len(dis) <= 1:
        return None

    dis_up = dis.upper()

    # Ignorar fuentes no aceptadas
    for src in IGNORE_SOURCES:
        if src in dis_up:
            return None

    # Extraer ID desnudo (quitar prefijo antes de : o _)
    naked = dis.split(':')[-1].split('_')[-1].strip().split(' ')[0]

    if 'OMIM' in dis_up:
        return f'OMIM:{naked}'
    elif 'ORPHANET' in dis_up or 'ORPH' in dis_up or 'ORDO' in dis_up or 'ORPHA' in dis_up:
        return f'ORPHA:{naked}'
    elif 'DOID' in dis_up or ('DO' in dis_up and 'ORDO' not in dis_up):
        return f'DOID:{naked}'
    elif 'HP' in dis_up or 'HPO' in dis_up:
        return f'HP:{naked}'
    elif 'EFO' in dis_up:
        return f'EFO:{naked}'
    elif 'UMLS' in dis_up or 'CUI' in dis_up or (len(dis) == 8 and dis.startswith('C')):
        return f'UMLS:{naked}'
    elif 'MESH' in dis_up or 'MSH' in dis_up or (dis[0].isalpha() and dis[1:].isnumeric()):
        return f'MESH:{naked}'
    elif dis.isnumeric():
        return f'MEDDRA:{naked}'
    else:
        return None


def process(raw_dir: str, out_dir: str, ctd_mapping: dict):
    os.makedirs(out_dir, exist_ok=True)
    infile = os.path.join(raw_dir, 'CTD_chemicals_diseases.tsv.gz')

    print(f'[ctdchemis] Procesando {infile} en streaming...')

    cau_rows = []
    trt_rows = []
    skipped  = 0

    with gzip.open(infile, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                skipped += 1
                continue

            chem_id    = cols[1].strip()   # CTD ChemicalID (ej: C046983)
            disease_id = cols[4].strip()   # DiseaseID (ej: MESH:D054198)
            evidence   = cols[5].strip().lower()

            # Solo evidencia directa
            if evidence not in CAUSES_EV and evidence not in TREATS_EV:
                skipped += 1
                continue

            # Mapear compuesto a InChIKey si esta disponible
            compound_id = ctd_mapping.get(chem_id, chem_id)

            # Normalizar ID de enfermedad
            norm_disease = normalize_disease_id(disease_id)
            if norm_disease is None:
                skipped += 1
                continue

            row = {
                'n1': compound_id,
                'n2': norm_disease,
            }

            if evidence in CAUSES_EV:
                cau_rows.append(row)
            else:
                trt_rows.append(row)

    cau = pd.DataFrame(cau_rows).drop_duplicates()
    trt = pd.DataFrame(trt_rows).drop_duplicates()

    cau_path = os.path.join(out_dir, 'CPD-cau-DIS.tsv')
    trt_path = os.path.join(out_dir, 'CPD-trt-DIS.tsv')

    cau.to_csv(cau_path, sep='\t', index=False)
    trt.to_csv(trt_path, sep='\t', index=False)

    print(f'[ctdchemis] CPD-cau-DIS: {len(cau):,} -> {cau_path}')
    print(f'[ctdchemis] CPD-trt-DIS: {len(trt):,} -> {trt_path}')
    print(f'[ctdchemis] Lineas descartadas: {skipped:,}')

    total_db = len(cau) + len(trt)
    stats_rows = []
    for rel, df in [('CPD-cau-DIS', cau), ('CPD-trt-DIS', trt)]:
        stats_rows.append({
            'database':       'ctdchemis',
            'relation':       rel,
            'gene_id_type':   'N/A',
            'assoc_id_type':  'MeSH Disease ID',
            'total_assoc':    len(df),
            'unique_genes':   df['n1'].nunique() if len(df) else 0,
            'unique_targets': df['n2'].nunique() if len(df) else 0,
            'db_total':       total_db,
        })

    pd.DataFrame(stats_rows).to_csv(
        os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False
    )
    print('[ctdchemis] Procesado completado.')


def main():
    args = parse_args()

    ctd_mapping = {}
    if args.ctd_mapping and os.path.exists(args.ctd_mapping):
        ctd_mapping = load_ctd_mapping(args.ctd_mapping)
    else:
        print('[ctdchemis] Sin mapeo CTD -> InChIKey: se usara ChemicalID de CTD directamente.')

    process(args.raw_dir, args.out_dir, ctd_mapping)
    sys.stderr.write('[ctdchemis] Done!\n')


if __name__ == '__main__':
    main()
