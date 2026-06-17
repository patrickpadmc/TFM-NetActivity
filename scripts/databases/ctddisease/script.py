#!/usr/bin/env python3
"""
ctddisease -- script.py
Procesado de CTD genes-diseases (CTD_genes_diseases.tsv.gz).

Logica (fiel a Bioteque):
  1. Saltar lineas sin DirectEvidence (inferidas)
  2. Mapear Entrez GeneID -> ENSG usando entrez2ensg.tsv
  3. Normalizar DiseaseID (estandarizar prefijo: MESH:, OMIM:, etc.)
  4. Deduplicar pares (ENSG, DiseaseID)

Columnas del archivo (linea 28):
  GeneSymbol | GeneID | DiseaseName | DiseaseID | DirectEvidence |
  InferenceChemicalName | InferenceScore | OmimIDs | PubMedIDs

Mapeo: Entrez -> ENSG via data/processed/entrez2ensg.tsv
  (ya construido en el proyecto, columnas: entrez_gene_id, ensembl_gene_id)

Output:
  processed/ctddisease/GEN-ass-DIS.tsv   (n1=ENSG, n2=DiseaseID)
  processed/ctddisease/stats.tsv

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_processed>
                   --entrez2ensg <ruta_entrez2ensg.tsv>
"""

import argparse
import gzip
import os
import sys
import pandas as pd


IGNORE_SOURCES = {'MEDGEN', 'SNOMED', 'BFO', 'ICD', 'NCI'}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir',      required=True)
    parser.add_argument('--out-dir',      required=True)
    parser.add_argument('--entrez2ensg',  required=True,
                        help='Ruta a data/processed/entrez2ensg.tsv')
    return parser.parse_args()


def load_entrez2ensg(path: str) -> dict:
    """Carga mapeo Entrez GeneID -> ENSG. Un Entrez puede mapear a varios ENSG."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    mapping = df.groupby('entrez_gene_id')['ensembl_gene_id'].apply(list).to_dict()
    print(f'[ctddisease] Mapeo Entrez->ENSG cargado: {len(mapping):,} entradas')
    return mapping


def normalize_disease_id(dis: str) -> str | None:
    """Estandariza prefijo del ID de enfermedad (replica parse_diseaseID de Bioteque)."""
    if not dis or len(dis) <= 1:
        return None

    dis_up = dis.upper()

    for src in IGNORE_SOURCES:
        if src in dis_up:
            return None

    naked = dis.split(':')[-1].split('_')[-1].strip().split(' ')[0]

    if 'OMIM' in dis_up:
        return f'OMIM:{naked}'
    elif 'ORPHANET' in dis_up or 'ORPH' in dis_up or 'ORDO' in dis_up or 'ORPHA' in dis_up:
        return f'ORPHA:{naked}'
    elif 'DOID' in dis_up:
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


def process(raw_dir: str, out_dir: str, entrez2ensg: dict):
    os.makedirs(out_dir, exist_ok=True)
    infile = os.path.join(raw_dir, 'CTD_genes_diseases.tsv.gz')

    print(f'[ctddisease] Procesando {infile} en streaming...')

    pairs = set()
    skipped_no_evidence = 0
    skipped_no_mapping  = 0
    skipped_no_disease  = 0
    n_lines             = 0

    with gzip.open(infile, 'rt', encoding='utf-8') as fh:
        # Saltar las 29 lineas de cabecera (igual que Bioteque)
        for _ in range(29):
            next(fh)

        for line in fh:
            n_lines += 1
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                skipped_no_evidence += 1
                continue

            entrez_id  = cols[1].strip()
            disease_id = cols[3].strip()
            evidence   = cols[4].strip()

            # Solo evidencia directa (igual que Bioteque)
            if not evidence:
                skipped_no_evidence += 1
                continue

            # Mapear Entrez -> ENSG
            ensg_list = entrez2ensg.get(entrez_id)
            if not ensg_list:
                skipped_no_mapping += 1
                continue

            # Normalizar ID de enfermedad
            norm_disease = normalize_disease_id(disease_id)
            if norm_disease is None:
                skipped_no_disease += 1
                continue

            for ensg in ensg_list:
                pairs.add((ensg, norm_disease))

    print(f'[ctddisease] Lineas procesadas:          {n_lines:,}')
    print(f'[ctddisease] Sin evidencia directa:      {skipped_no_evidence:,}')
    print(f'[ctddisease] Sin mapeo Entrez->ENSG:     {skipped_no_mapping:,}')
    print(f'[ctddisease] DiseaseID descartado:       {skipped_no_disease:,}')
    print(f'[ctddisease] Pares unicos (ENSG, DIS):   {len(pairs):,}')

    # Guardar
    out_path = os.path.join(out_dir, 'GEN-ass-DIS.tsv')
    with open(out_path, 'w') as f:
        f.write('n1\tn2\n')
        for ensg, dis in sorted(pairs):
            f.write(f'{ensg}\t{dis}\n')

    print(f'[ctddisease] GEN-ass-DIS -> {out_path}')

    # Stats
    df = pd.read_csv(out_path, sep='\t')
    stats = pd.DataFrame([{
        'database':       'ctddisease',
        'relation':       'GEN-ass-DIS',
        'gene_id_type':   'ENSG',
        'assoc_id_type':  'MeSH/OMIM/ORPHA Disease ID',
        'total_assoc':    len(df),
        'unique_genes':   df['n1'].nunique(),
        'unique_targets': df['n2'].nunique(),
        'db_total':       len(df),
    }])
    stats.to_csv(os.path.join(out_dir, 'stats.tsv'), sep='\t', index=False)
    print('[ctddisease] Procesado completado.')


def main():
    args = parse_args()
    entrez2ensg = load_entrez2ensg(args.entrez2ensg)
    process(args.raw_dir, args.out_dir, entrez2ensg)
    sys.stderr.write('[ctddisease] Done!\n')


if __name__ == '__main__':
    main()
