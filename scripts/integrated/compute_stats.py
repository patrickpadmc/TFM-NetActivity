#!/usr/bin/env python3
"""
compute_stats.py
Genera 6 tablas de estadisticas para todas las databases procesadas.
Output: TSV por tabla + un HTML interactivo con las 6 tablas.

Definiciones:
  SF (Sin Filtrar): antes del filtro de calidad, leyendo archivos raw
  F  (Filtrado):    despues del filtro de calidad, leyendo archivos procesados
  OG (Original):    IDs tal como estan en los archivos raw/procesados
  HOM (Homogenizado): genes -> ENSG; asoc -> DOID (enf) / BTO (tej) / InChIKey (cpd)

Tablas generadas:
  1. IDs originales, SF y F
  2. IDs homogenizados, SF y F
  3. Comparacion SF: original vs homogenizado
  4. Comparacion F: original vs homogenizado
  5. Resumen filtrado + homogenizado
  6. Desglose por relacion

Uso:
  python compute_stats.py [--base-dir PATH] [--out-dir PATH]
"""

import argparse, gzip, glob, os, sys, math
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Dict, Set
from collections import defaultdict

# ─── Rutas por defecto ────────────────────────────────────────────────────────
BASE_DIR     = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity"
BIOTEQUE_MAP = "/beegfs/home/ppadmoremcc/work/external/bioteque/metadata/mappings"

# ─── Parametros de calidad (identicos a cada script) ─────────────────────────
OT_SIGN      = 0.7
OT_MIN_SCORE = 0.3
OT_MIN_HITS  = 100
STR_MIN      = 700
HPA_RNA_MAX  = 250
HPA_RNA_POS  = 0.5
HPA_RNA_NEG  = -0.5
DOROTHEA_AB  = {'A', 'B'}
DOROTHEA_CD  = {'C', 'D'}
CTD_EV       = {'marker/mechanism', 'marker', 'mechanism', 'therapeutic'}
COEX_PCT     = 0.1

# ─── Estructura de datos ──────────────────────────────────────────────────────
@dataclass
class Counts:
    total:     Optional[int] = None
    gene_unq:  Optional[int] = None
    assoc_unq: Optional[int] = None

@dataclass
class RelStats:
    db:          str = ''
    relation:    str = ''
    gene_id_og:  str = ''
    assoc_id_og: str = ''
    sf_og:  Counts = field(default_factory=Counts)
    sf_hom: Counts = field(default_factory=Counts)
    f_og:   Counts = field(default_factory=Counts)
    f_hom:  Counts = field(default_factory=Counts)


def cnt(pairs_set: set) -> Counts:
    """Convierte un set de (gene, assoc) a Counts."""
    if pairs_set is None:
        return Counts()
    total = len(pairs_set)
    genes = len({p[0] for p in pairs_set})
    assoc = len({p[1] for p in pairs_set})
    return Counts(total, genes, assoc)


def cnt_from_df(df: Optional[pd.DataFrame]) -> Counts:
    if df is None or len(df) == 0:
        return Counts()
    return Counts(len(df), df['n1'].nunique(), df['n2'].nunique())


def read_proc(proc_dir: str, db: str, rel: str) -> Optional[pd.DataFrame]:
    path = os.path.join(proc_dir, db, f"{rel}.tsv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep='\t', low_memory=False)
    if 'n1' not in df.columns or 'n2' not in df.columns:
        return None
    return df


# ─── Carga de mappings ────────────────────────────────────────────────────────

def load_mappings(base_dir: str) -> Dict:
    maps = {}

    # Entrez -> ENSG
    p = f"{base_dir}/data/processed/entrez2ensg.tsv"
    if os.path.exists(p):
        df = pd.read_csv(p, sep='\t', dtype=str)
        maps['entrez2ensg'] = dict(zip(df['entrez_gene_id'], df['ensembl_gene_id']))
        print(f"[maps] entrez2ensg: {len(maps['entrez2ensg']):,}")
    else:
        maps['entrez2ensg'] = {}

    # ENSP -> ENSG y UniProt -> ENSG
    p = f"{base_dir}/data/processed/ensp2ensg.tsv.gz"
    ensp2ensg = {}
    uniprot2ensg = {}
    if os.path.exists(p):
        with gzip.open(p, 'rt') as f:
            next(f)
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 5: continue
                ensg, ensp, uniprot, db_name = cols[0], cols[2], cols[3], cols[4]
                if ensg.startswith('ENSG') and ensp.startswith('ENSP'):
                    if ensp not in ensp2ensg:
                        ensp2ensg[ensp] = ensg
                if db_name == 'Uniprot/SWISSPROT' and ensg.startswith('ENSG') and uniprot:
                    if uniprot not in uniprot2ensg:
                        uniprot2ensg[uniprot] = ensg
    maps['ensp2ensg']    = ensp2ensg
    maps['uniprot2ensg'] = uniprot2ensg
    print(f"[maps] ensp2ensg: {len(ensp2ensg):,} | uniprot2ensg: {len(uniprot2ensg):,}")

    # Symbol -> ENSG (HGNC)
    p = f"{base_dir}/data/processed/symbol2ensg.tsv"
    if os.path.exists(p):
        df = pd.read_csv(p, sep='\t', low_memory=False, usecols=['symbol', 'ensembl_gene_id'])
        df = df.dropna(subset=['ensembl_gene_id'])
        maps['symbol2ensg'] = dict(zip(df['symbol'], df['ensembl_gene_id']))
        print(f"[maps] symbol2ensg: {len(maps['symbol2ensg']):,}")
    else:
        maps['symbol2ensg'] = {}

    # Tejidos -> BTO (Bioteque)
    p = f"{BIOTEQUE_MAP}/TIS/tis2id.tsv"
    tis2bto = {}
    if os.path.exists(p):
        with open(p, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    tis2bto[parts[0].strip().upper()] = parts[1].strip()
    maps['tis2bto'] = tis2bto
    print(f"[maps] tis2bto: {len(tis2bto):,}")

    # Enfermedades -> DOID (Bioteque doid.tsv)
    p = f"{BIOTEQUE_MAP}/DIS/doid.tsv"
    dis2doid = {}
    if os.path.exists(p):
        with open(p, 'r') as f:
            next(f)
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[0] not in dis2doid:
                    dis2doid[parts[0].strip()] = parts[1].strip()
    maps['dis2doid'] = dis2doid
    print(f"[maps] dis2doid: {len(dis2doid):,}")

    # CTD ID -> InChIKey
    p = f"{BIOTEQUE_MAP}/CPD/ctd.tsv"
    ctd2ikey = {}
    if os.path.exists(p):
        with open(p) as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 3 and cols[2]:
                    ctd2ikey[cols[0]] = cols[2]
    maps['ctd2ikey'] = ctd2ikey
    print(f"[maps] ctd2ikey: {len(ctd2ikey):,}")

    # Repohub pert_iname -> InChIKey
    p = f"{BIOTEQUE_MAP}/CPD/repohub.tsv"
    repo2ikey = {}
    if os.path.exists(p):
        with open(p) as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 2 and cols[1]:
                    repo2ikey[cols[0].strip('"')] = cols[1]
    maps['repo2ikey'] = repo2ikey
    print(f"[maps] repo2ikey: {len(repo2ikey):,}")

    return maps


# ─── Normalizacion de IDs ─────────────────────────────────────────────────────
IGNORE_DIS = {'MEDGEN', 'SNOMED', 'BFO', 'ICD', 'NCI'}

def norm_disease(dis: str) -> Optional[str]:
    if not dis or len(dis) <= 1: return None
    u = dis.upper()
    for s in IGNORE_DIS:
        if s in u: return None
    naked = dis.split(':')[-1].split('_')[-1].strip().split(' ')[0]
    if 'OMIM' in u:    return f'OMIM:{naked}'
    if 'ORPH' in u:    return f'ORPHA:{naked}'
    if 'DOID' in u:    return f'DOID:{naked}'
    if 'HP' in u:      return f'HP:{naked}'
    if 'EFO' in u:     return f'EFO:{naked}'
    if 'UMLS' in u or (len(dis)==8 and dis.startswith('C')): return f'UMLS:{naked}'
    if 'MESH' in u or (dis[0].isalpha() and dis[1:].isnumeric()): return f'MESH:{naked}'
    return None

def hom_disease(dis_id: str, dis2doid: dict) -> Optional[str]:
    if not dis_id: return None
    return dis2doid.get(dis_id, dis_id)

def hom_tissue(tis_name: str, tis2bto: dict) -> Optional[str]:
    return tis2bto.get(tis_name.upper(), None)

def norm_ot_disease(dis_id: str) -> str:
    if dis_id.startswith('Orphanet_'):
        return 'ORPHA:' + dis_id.split('_', 1)[1]
    return dis_id.replace('_', ':', 1)


# ─── Stats por database ───────────────────────────────────────────────────────

def stats_opentargets(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_dir = f"{base_dir}/data/raw/databases/opentargets"
    parquets = sorted(glob.glob(f"{raw_dir}/*.parquet"))
    r = RelStats(db='opentargets', relation='GEN-ass-DIS',
                 gene_id_og='ENSG', assoc_id_og='EFO/DOID/ORPHA/HP')

    if parquets:
        print("[stats] opentargets SF...")
        sf_og, sf_hom = set(), set()
        for p in parquets:
            df = pd.read_parquet(p, columns=['diseaseId', 'targetId'])
            for ensg, dis in zip(df['targetId'], df['diseaseId']):
                dis_og  = norm_ot_disease(dis)
                dis_hom = hom_disease(dis_og, maps['dis2doid'])
                sf_og.add((ensg, dis_og))
                if dis_hom: sf_hom.add((ensg, dis_hom))
        r.sf_og  = cnt(sf_og)
        r.sf_hom = cnt(sf_hom)

    # F from processed
    df_f = read_proc(proc_dir, 'opentargets', 'GEN-ass-DIS')
    r.f_og = cnt_from_df(df_f)
    if df_f is not None:
        f_hom = set()
        for ensg, dis in zip(df_f['n1'], df_f['n2']):
            dis_hom = hom_disease(dis, maps['dis2doid'])
            if dis_hom: f_hom.add((ensg, dis_hom))
        r.f_hom = cnt(f_hom)
    return [r]


def stats_hpa_proteome(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/hpa_proteome/normal_ihc_data.tsv"
    results = []

    for rel, level_filter in [('GEN-pab-TIS', {'High'}), ('GEN-pdf-TIS', {'Low'})]:
        r = RelStats(db='hpa_proteome', relation=rel,
                     gene_id_og='ENSG', assoc_id_og='Tissue name')
        if os.path.exists(raw_path):
            print(f"[stats] hpa_proteome {rel} SF...")
            df = pd.read_csv(raw_path, sep='\t', low_memory=False,
                             usecols=['Gene', 'Tissue', 'Level', 'Reliability'])
            # SF: all rows excluding Uncertain (any level)
            df_sf = df[df['Reliability'] != 'Uncertain']
            sf_og  = set(zip(df_sf['Gene'], df_sf['Tissue']))
            sf_hom_set = set()
            for ensg, tis in sf_og:
                bto = hom_tissue(tis, maps['tis2bto'])
                if bto: sf_hom_set.add((ensg, bto))
            r.sf_og  = cnt(sf_og)
            r.sf_hom = cnt(sf_hom_set)

        df_f = read_proc(proc_dir, 'hpa_proteome', rel)
        r.f_og = cnt_from_df(df_f)
        if df_f is not None:
            f_hom = set()
            for ensg, tis in zip(df_f['n1'], df_f['n2']):
                bto = hom_tissue(tis, maps['tis2bto'])
                if bto: f_hom.add((ensg, bto))
            r.f_hom = cnt(f_hom)
        results.append(r)
    return results


def stats_hpa_rna_cons(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/hpa_rna_cons/rna_tissue_consensus.tsv"
    results = []
    for rel in ['GEN-upr-TIS', 'GEN-dwr-TIS']:
        r = RelStats(db='hpa_rna_cons', relation=rel,
                     gene_id_og='ENSG', assoc_id_og='Tissue name')
        if os.path.exists(raw_path):
            print(f"[stats] hpa_rna_cons {rel} SF...")
            df = pd.read_csv(raw_path, sep='\t', low_memory=False,
                             usecols=['Gene', 'Tissue'])
            sf_og  = set(zip(df['Gene'], df['Tissue']))
            sf_hom = set()
            for ensg, tis in sf_og:
                bto = hom_tissue(tis, maps['tis2bto'])
                if bto: sf_hom.add((ensg, bto))
            r.sf_og  = cnt(sf_og)
            r.sf_hom = cnt(sf_hom)

        df_f = read_proc(proc_dir, 'hpa_rna_cons', rel)
        r.f_og = cnt_from_df(df_f)
        if df_f is not None:
            f_hom = set()
            for ensg, tis in zip(df_f['n1'], df_f['n2']):
                bto = hom_tissue(tis, maps['tis2bto'])
                if bto: f_hom.add((ensg, bto))
            r.f_hom = cnt(f_hom)
        results.append(r)
    return results


def stats_ctddisease(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/ctddisease/CTD_genes_diseases.tsv.gz"
    r = RelStats(db='ctddisease', relation='GEN-ass-DIS',
                 gene_id_og='Entrez', assoc_id_og='MeSH/OMIM')
    if os.path.exists(raw_path):
        print("[stats] ctddisease SF (streaming)...")
        sf_genes, sf_dis = set(), set()
        f_genes_og, f_dis_og = set(), set()
        f_genes_hom, f_dis_hom = set(), set()
        sf_total = 0
        f_total_og = 0
        e2e = maps['entrez2ensg']
        with gzip.open(raw_path, 'rt', encoding='utf-8') as fh:
            for _ in range(29): next(fh)
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 5: continue
                entrez = cols[1].strip()
                dis_raw = cols[3].strip()
                evidence = cols[4].strip()
                dis_og = norm_disease(dis_raw)
                if dis_og is None: continue
                # SF: all rows con cualquier evidence
                sf_genes.add(entrez)
                sf_dis.add(dis_og)
                sf_total += 1
                # F: solo DirectEvidence
                if evidence:
                    ensg = e2e.get(entrez)
                    f_genes_og.add(entrez)
                    f_dis_og.add(dis_og)
                    f_total_og += 1
                    if ensg:
                        dis_hom = hom_disease(dis_og, maps['dis2doid'])
                        f_genes_hom.add(ensg)
                        if dis_hom: f_dis_hom.add(dis_hom)

        r.sf_og = Counts(sf_total, len(sf_genes), len(sf_dis))
        # SF HOM: genes mapeados a ENSG
        sf_genes_hom = {e2e[g] for g in sf_genes if g in e2e}
        sf_dis_hom   = {hom_disease(d, maps['dis2doid']) for d in sf_dis}
        sf_dis_hom.discard(None)
        r.sf_hom = Counts(sf_total, len(sf_genes_hom), len(sf_dis_hom))
        r.f_og   = Counts(f_total_og, len(f_genes_og), len(f_dis_og))
        r.f_hom  = Counts(f_total_og, len(f_genes_hom), len(f_dis_hom))
    else:
        df_f = read_proc(proc_dir, 'ctddisease', 'GEN-ass-DIS')
        r.f_og = cnt_from_df(df_f)
    return [r]


def stats_ctdchemis(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/ctdchemis/CTD_chemicals_diseases.tsv.gz"
    results = []
    ctd2ikey = maps['ctd2ikey']

    for rel, ev_set in [('CPD-cau-DIS', {'marker/mechanism','marker','mechanism'}),
                        ('CPD-trt-DIS', {'therapeutic'})]:
        r = RelStats(db='ctdchemis', relation=rel,
                     gene_id_og='N/A', assoc_id_og='MeSH Disease')
        if os.path.exists(raw_path):
            print(f"[stats] ctdchemis {rel} SF (streaming)...")
            sf_cpd, sf_dis = set(), set()
            f_cpd_og, f_dis_og = set(), set()
            f_cpd_hom, f_dis_hom = set(), set()
            sf_total = f_total = 0
            with gzip.open(raw_path, 'rt', encoding='utf-8') as fh:
                for line in fh:
                    if line.startswith('#'): continue
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 6: continue
                    chem_id = cols[1].strip()
                    dis_raw = cols[4].strip()
                    evidence = cols[5].strip().lower()
                    dis_og = norm_disease(dis_raw)
                    if dis_og is None: continue
                    # SF: todas las DirectEvidence (cualquier tipo)
                    if evidence in CTD_EV:
                        sf_cpd.add(chem_id)
                        sf_dis.add(dis_og)
                        sf_total += 1
                    # F: solo las del relation actual
                    if evidence in ev_set:
                        f_cpd_og.add(chem_id)
                        f_dis_og.add(dis_og)
                        f_total += 1
                        ikey = ctd2ikey.get(chem_id, chem_id)
                        f_cpd_hom.add(ikey)
                        dis_hom = hom_disease(dis_og, maps['dis2doid'])
                        if dis_hom: f_dis_hom.add(dis_hom)

            sf_cpd_hom = {ctd2ikey.get(c, c) for c in sf_cpd}
            sf_dis_hom = {hom_disease(d, maps['dis2doid']) for d in sf_dis}
            sf_dis_hom.discard(None)
            r.sf_og  = Counts(sf_total, len(sf_cpd), len(sf_dis))
            r.sf_hom = Counts(sf_total, len(sf_cpd_hom), len(sf_dis_hom))
            r.f_og   = Counts(f_total, len(f_cpd_og), len(f_dis_og))
            r.f_hom  = Counts(f_total, len(f_cpd_hom), len(f_dis_hom))
        else:
            df_f = read_proc(proc_dir, 'ctdchemis', rel)
            r.f_og = cnt_from_df(df_f)
        results.append(r)
    return results


def stats_jensentissuecurated(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/jensentissuecurated/human_tissue_knowledge_full.tsv"
    r = RelStats(db='jensentissuecurated', relation='GEN-ass-TIS',
                 gene_id_og='ENSP', assoc_id_og='BTO ID')
    if os.path.exists(raw_path):
        print("[stats] jensentissuecurated SF...")
        sf_og, sf_hom = set(), set()
        ensp2e = maps['ensp2ensg']
        with open(raw_path, 'r') as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 3: continue
                ensp   = cols[0].strip()
                bto_id = cols[2].strip()
                if bto_id == 'BTO:0000000': continue
                sf_og.add((ensp, bto_id))
                ensg = ensp2e.get(ensp)
                if ensg: sf_hom.add((ensg, bto_id))
        r.sf_og  = cnt(sf_og)
        r.sf_hom = cnt(sf_hom)

    df_f = read_proc(proc_dir, 'jensentissuecurated', 'GEN-ass-TIS')
    r.f_og  = cnt_from_df(df_f)
    r.f_hom = cnt_from_df(df_f)  # ya tiene ENSG y BTO
    return [r]


def stats_reactome(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/reactome/UniProt2Reactome_All_Levels.txt"
    r = RelStats(db='reactome', relation='GEN-ass-PWY',
                 gene_id_og='UniProt', assoc_id_og='Reactome R-HSA-')
    if os.path.exists(raw_path):
        print("[stats] reactome SF...")
        sf_all_og, sf_hsa_og = set(), set()
        sf_hsa_hom = set()
        up2e = maps['uniprot2ensg']
        with open(raw_path, 'r') as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 2: continue
                uniprot = cols[0].strip()
                pwy     = cols[1].strip()
                sf_all_og.add((uniprot, pwy))
                if pwy.startswith('R-HSA-'):
                    sf_hsa_og.add((uniprot, pwy))
                    ensg = up2e.get(uniprot)
                    if ensg: sf_hsa_hom.add((ensg, pwy))
        # SF = all rows (any species); F = only R-HSA- (human)
        r.sf_og  = cnt(sf_all_og)
        r.sf_hom = cnt(sf_hsa_hom)  # SF HOM: only human + ENSG
        r.f_og   = cnt(sf_hsa_og)
        r.f_hom  = cnt(sf_hsa_hom)
    else:
        df_f = read_proc(proc_dir, 'reactome', 'GEN-ass-PWY')
        r.f_og = r.f_hom = cnt_from_df(df_f)
    return [r]


def stats_repohub(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/repohub/repo-drug-annotation-20200324.txt"
    r = RelStats(db='repohub', relation='CPD-int-GEN',
                 gene_id_og='Gene symbol', assoc_id_og='pert_iname')
    s2e = maps['symbol2ensg']
    r2k = maps['repo2ikey']

    if os.path.exists(raw_path):
        print("[stats] repohub SF...")
        sf_og, sf_hom = set(), set()
        header = None
        with open(raw_path, 'r') as f:
            for line in f:
                if line.startswith('!'): continue
                cols = line.rstrip('\n').split('\t')
                if header is None: header = cols; continue
                d = dict(zip(header, cols))
                cpd = d.get('pert_iname', '').strip()
                tgt = d.get('target', '').strip()
                if not tgt: continue
                for sym in tgt.split('|'):
                    sym = sym.strip()
                    if not sym: continue
                    sf_og.add((sym, cpd))
                    ensg = s2e.get(sym)
                    ikey = r2k.get(cpd, cpd)
                    if ensg: sf_hom.add((ensg, ikey))
        r.sf_og  = cnt(sf_og)
        r.sf_hom = cnt(sf_hom)

    df_f = read_proc(proc_dir, 'repohub', 'CPD-int-GEN')
    r.f_og = cnt_from_df(df_f)
    if df_f is not None:
        f_hom = set()
        for cpd, ensg in zip(df_f['n1'], df_f['n2']):
            ikey = r2k.get(cpd, cpd)
            f_hom.add((ensg, ikey))
        r.f_hom = cnt(f_hom)
    return [r]


def stats_string(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    raw_path = f"{base_dir}/data/raw/databases/string/9606.protein.links.full.v12.0.txt.gz"
    r = RelStats(db='string', relation='GEN-ppi-GEN',
                 gene_id_og='ENSP', assoc_id_og='ENSP')
    if os.path.exists(raw_path):
        print("[stats] string SF (streaming)...")
        sf_genes, f_genes = set(), set()
        sf_total = f_total = 0
        e2e = maps['ensp2ensg']
        with gzip.open(raw_path, 'rt') as fh:
            fh.readline()
            for line in fh:
                cols = line.rstrip('\n').split(' ')
                if len(cols) < 3: continue
                score = int(cols[-1])
                ensp1 = cols[0].split('.')[-1]
                ensp2 = cols[1].split('.')[-1]
                sf_genes.update([ensp1, ensp2])
                sf_total += 1
                if score >= STR_MIN:
                    ensg1 = e2e.get(ensp1)
                    ensg2 = e2e.get(ensp2)
                    if ensg1: f_genes.add(ensg1)
                    if ensg2: f_genes.add(ensg2)
                    f_total += 1
        sf_genes_hom = {e2e[g] for g in sf_genes if g in e2e}
        r.sf_og  = Counts(sf_total, len(sf_genes), len(sf_genes))
        r.sf_hom = Counts(sf_total, len(sf_genes_hom), len(sf_genes_hom))
        r.f_og   = Counts(f_total, len(f_genes), len(f_genes))
        r.f_hom  = Counts(f_total, len(f_genes), len(f_genes))
    else:
        df_f = read_proc(proc_dir, 'string', 'GEN-ppi-GEN')
        r.f_og = r.f_hom = cnt_from_df(df_f)
    return [r]


def stats_omnipath(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    int_path = f"{base_dir}/data/raw/databases/omnipath/omnipath_webservice_interactions__latest.tsv.gz"
    enz_path = f"{base_dir}/data/raw/databases/omnipath/omnipath_webservice_enz_sub__latest.tsv.gz"
    up2e = maps['uniprot2ensg']
    results = []

    # PPI
    r_ppi = RelStats(db='omnipath', relation='GEN-ppi-GEN',
                     gene_id_og='UniProt', assoc_id_og='UniProt')
    if os.path.exists(int_path):
        print("[stats] omnipath PPI SF...")
        df = pd.read_csv(int_path, sep='\t', low_memory=False)
        reviewed = set(up2e.keys())
        # SF: all human rows
        sf_df = df[(df['ncbi_tax_id_source']==9606) & (df['ncbi_tax_id_target']==9606)]
        sf_pairs = set(zip(sf_df['source'], sf_df['target']))
        r_ppi.sf_og = cnt(sf_pairs)
        sf_hom = set()
        for s, t in sf_pairs:
            e1, e2 = up2e.get(s), up2e.get(t)
            if e1 and e2 and e1 != e2: sf_hom.add(tuple(sorted([e1,e2])))
        r_ppi.sf_hom = cnt(sf_hom)
    df_f = read_proc(proc_dir, 'omnipath', 'GEN-ppi-GEN')
    r_ppi.f_og = r_ppi.f_hom = cnt_from_df(df_f)
    results.append(r_ppi)

    # Phosphorylation + dephosphorylation from enz_sub
    if os.path.exists(enz_path):
        print("[stats] omnipath PTM SF...")
        df = pd.read_csv(enz_path, sep='\t', low_memory=False)
        if 'ncbi_tax_id' in df.columns:
            df_h = df[df['ncbi_tax_id'] == 9606]
        else:
            df_h = df
        for rel, mod in [('GEN-pho-GEN', 'phosphorylation'),
                         ('GEN-dph-GEN', 'dephosphorylation')]:
            r = RelStats(db='omnipath', relation=rel,
                         gene_id_og='UniProt', assoc_id_og='UniProt')
            sf_df = df_h
            sf_mod = sf_df[sf_df['modification'] == mod]
            sf_pairs = set(zip(sf_mod['enzyme'], sf_mod['substrate']))
            r.sf_og = cnt(sf_pairs)
            sf_hom = set()
            for s, t in sf_pairs:
                e1, e2 = up2e.get(s), up2e.get(t)
                if e1 and e2 and e1 != e2: sf_hom.add((e1, e2))
            r.sf_hom = cnt(sf_hom)
            df_f = read_proc(proc_dir, 'omnipath', rel)
            r.f_og = r.f_hom = cnt_from_df(df_f)
            results.append(r)
    return results


def stats_dorothea(base_dir: str, proc_dir: str, maps: dict,
                   levels: set, db_name: str) -> List[RelStats]:
    int_path = f"{base_dir}/data/raw/databases/omnipath/omnipath_webservice_interactions__latest.tsv.gz"
    up2e = maps['uniprot2ensg']
    results = []

    df = None
    if os.path.exists(int_path):
        print(f"[stats] {db_name} SF...")
        df_all = pd.read_csv(int_path, sep='\t', low_memory=False)
        df_dor = df_all[df_all['dorothea'] == True]
        df_h   = df_dor[(df_dor['ncbi_tax_id_source']==9606) & (df_dor['ncbi_tax_id_target']==9606)]
        # Level filter
        level_vals = set()
        for x in df_h['dorothea_level'].dropna().unique():
            for lvl in levels:
                if lvl in x: level_vals.add(x)
        df = df_h[df_h['dorothea_level'].isin(level_vals)]

    for rel, col in [('GEN-reg-GEN', None),
                     ('GEN-upr-GEN', 'consensus_stimulation'),
                     ('GEN-dwr-GEN', 'consensus_inhibition')]:
        r = RelStats(db=db_name, relation=rel,
                     gene_id_og='UniProt', assoc_id_og='UniProt')
        if df is not None:
            sub = df if col is None else df[df[col] == True]
            sf_pairs = set(zip(sub['source'], sub['target']))
            r.sf_og = cnt(sf_pairs)
            sf_hom = set()
            for s, t in sf_pairs:
                e1, e2 = up2e.get(s), up2e.get(t)
                if e1 and e2 and e1 != e2: sf_hom.add((e1, e2))
            r.sf_hom = cnt(sf_hom)
        df_f = read_proc(proc_dir, db_name, rel)
        r.f_og = r.f_hom = cnt_from_df(df_f)
        results.append(r)
    return results


def stats_coexpressdb(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    r = RelStats(db='coexpressdb', relation='GEN-cex-GEN',
                 gene_id_og='Entrez', assoc_id_og='Entrez')
    # SF: count all lines across coex_data files (expensive but correct)
    coex_dir = f"{base_dir}/data/raw/databases/coexpressdb/coex_data"
    if os.path.isdir(coex_dir):
        print("[stats] coexpressdb SF (counting all lines)...")
        all_genes_files = os.listdir(coex_dir)
        n_files = len(all_genes_files)
        cutoff = math.ceil(n_files * COEX_PCT)
        sf_total = 0
        sf_source = set(all_genes_files)
        sf_target = set()
        # Only count source genes, and estimate targets
        for fname in all_genes_files:
            fpath = os.path.join(coex_dir, fname)
            with open(fpath, 'r') as f:
                for line in f:
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 2: continue
                    sf_total += 1
                    sf_target.add(cols[0].strip())
        sf_all_genes = sf_source | sf_target
        e2e = maps['entrez2ensg']
        sf_genes_hom = {e2e[g] for g in sf_all_genes if g in e2e}
        r.sf_og  = Counts(sf_total, len(sf_source), len(sf_target))
        r.sf_hom = Counts(sf_total, len(sf_genes_hom), len(sf_genes_hom))

    df_f = read_proc(proc_dir, 'coexpressdb', 'GEN-cex-GEN')
    r.f_og = r.f_hom = cnt_from_df(df_f)
    return [r]


def stats_ccle_family(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    """cclerna, cclemut, cclecnv."""
    results = []
    s2e = maps['symbol2ensg']

    # cclerna SF
    expr_path = f"{base_dir}/data/raw/databases/cclerna_HMZ/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    model_path = f"{base_dir}/data/raw/databases/cclerna_HMZ/Model.csv"
    for rel in ['CLL-upr-GEN', 'CLL-dwr-GEN']:
        r = RelStats(db='cclerna_HMZ', relation=rel,
                     gene_id_og='Gene symbol', assoc_id_og='CVCL')
        if os.path.exists(expr_path) and os.path.exists(model_path):
            print(f"[stats] cclerna_HMZ {rel} SF...")
            model_df = pd.read_csv(model_path, usecols=['ModelID','RRID'], low_memory=False)
            model_df = model_df.dropna(subset=['RRID'])
            model_df = model_df[model_df['RRID'].str.startswith('CVCL_')]
            m2cvcl = dict(zip(model_df['ModelID'], model_df['RRID']))
            df = pd.read_csv(expr_path, index_col=0, low_memory=False, nrows=0)
            cols = df.columns.tolist()
            n_genes = len(cols)
            n_cls   = sum(1 for x in open(expr_path) if x.strip()) - 1
            ensg_set = {s2e[c.split(' ')[0]] for c in cols if c.split(' ')[0] in s2e}
            cvcl_set = set(m2cvcl.values())
            sf_total = n_genes * min(n_cls, len(m2cvcl))
            r.sf_og  = Counts(sf_total, len(cvcl_set), n_genes)
            r.sf_hom = Counts(sf_total, len(cvcl_set), len(ensg_set))
        df_f = read_proc(proc_dir, 'cclerna_HMZ', rel)
        r.f_og = r.f_hom = cnt_from_df(df_f)
        results.append(r)

    # cclemut SF: already binary, SF ~ F
    r = RelStats(db='cclemut_HMZ', relation='CLL-mut-GEN',
                 gene_id_og='Gene symbol', assoc_id_og='CVCL')
    df_f = read_proc(proc_dir, 'cclemut_HMZ', 'CLL-mut-GEN')
    r.f_og = r.f_hom = r.sf_og = r.sf_hom = cnt_from_df(df_f)
    results.append(r)

    # cclecnv
    cnv_path = f"{base_dir}/data/raw/databases/cclecnv_HMZ/gene_attribute_edges.txt.gz"
    for rel, w_filter in [('CLL-cnu-GEN', 1.0), ('CLL-cnd-GEN', -1.0)]:
        r = RelStats(db='cclecnv_HMZ', relation=rel,
                     gene_id_og='Gene symbol', assoc_id_og='Cell line name')
        if os.path.exists(cnv_path):
            print(f"[stats] cclecnv_HMZ {rel} SF...")
            sf_og, sf_hom = set(), set()
            with gzip.open(cnv_path, 'rt') as f:
                next(f); next(f)
                for line in f:
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 7: continue
                    sym = cols[0].strip()
                    cl  = cols[3].strip()
                    try:
                        w = float(cols[6])
                    except ValueError:
                        continue
                    if w == w_filter:
                        sf_og.add((cl, sym))
                        ensg = s2e.get(sym)
                        if ensg: sf_hom.add((cl, ensg))
            r.sf_og  = cnt(sf_og)
            r.sf_hom = cnt(sf_hom)
        df_f = read_proc(proc_dir, 'cclecnv_HMZ', rel)
        r.f_og = r.f_hom = cnt_from_df(df_f)
        results.append(r)

    return results


def stats_achilles(base_dir: str, proc_dir: str, maps: dict) -> List[RelStats]:
    ach_path = f"{base_dir}/data/raw/databases/achilles_HMZ/gene_attribute_edges.txt.gz"
    s2e = maps['symbol2ensg']
    results = []
    for rel, w_filter in [('CLL-bfn-GEN', -1.0), ('CLL-gfn-GEN', 1.0)]:
        r = RelStats(db='achilles_HMZ', relation=rel,
                     gene_id_og='Gene symbol', assoc_id_og='Cell line name')
        if os.path.exists(ach_path):
            print(f"[stats] achilles_HMZ {rel} SF...")
            sf_og, sf_hom = set(), set()
            with gzip.open(ach_path, 'rt') as f:
                next(f); next(f)
                for line in f:
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) < 7: continue
                    sym = cols[0].strip()
                    cl  = cols[3].strip()
                    try:
                        w = float(cols[6])
                    except ValueError:
                        continue
                    if w == w_filter:
                        sf_og.add((cl, sym))
                        ensg = s2e.get(sym)
                        if ensg: sf_hom.add((cl, ensg))
            r.sf_og  = cnt(sf_og)
            r.sf_hom = cnt(sf_hom)
        df_f = read_proc(proc_dir, 'achilles_HMZ', rel)
        r.f_og = r.f_hom = cnt_from_df(df_f)
        results.append(r)
    return results


# ─── Construccion de tablas ───────────────────────────────────────────────────

def fmt(v: Optional[int]) -> str:
    if v is None: return '~'
    return f'{v:,}'

def dif(a: Optional[int], b: Optional[int]) -> str:
    if a is None or b is None: return '~'
    return f'{a-b:,}'


def build_table1(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 1: IDs originales, SF y F."""
    data = []
    for r in rows:
        data.append({
            'Database':         r.db,
            'Relation':         r.relation,
            'Gene ID OG':       r.gene_id_og,
            'Assoc ID OG':      r.assoc_id_og,
            'ASOC TOT OG SF':   fmt(r.sf_og.total),
            'GENE UNQ OG SF':   fmt(r.sf_og.gene_unq),
            'ASOC UNQ OG SF':   fmt(r.sf_og.assoc_unq),
            'ASOC TOT OG F':    fmt(r.f_og.total),
            'GENE UNQ OG F':    fmt(r.f_og.gene_unq),
            'ASOC UNQ OG F':    fmt(r.f_og.assoc_unq),
            'DIF ASOC TOT':     dif(r.sf_og.total,     r.f_og.total),
            'DIF GENE UNQ':     dif(r.sf_og.gene_unq,  r.f_og.gene_unq),
            'DIF ASOC UNQ':     dif(r.sf_og.assoc_unq, r.f_og.assoc_unq),
        })
    return pd.DataFrame(data)


def build_table2(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 2: IDs homogenizados (ENSG + vocab estandar), SF y F."""
    data = []
    for r in rows:
        data.append({
            'Database':          r.db,
            'Relation':          r.relation,
            'Gene ID HOM':       'ENSG',
            'Assoc ID HOM':      'DOID/BTO/InChIKey',
            'ASOC TOT HOM SF':   fmt(r.sf_hom.total),
            'GENE UNQ HOM SF':   fmt(r.sf_hom.gene_unq),
            'ASOC UNQ HOM SF':   fmt(r.sf_hom.assoc_unq),
            'ASOC TOT HOM F':    fmt(r.f_hom.total),
            'GENE UNQ HOM F':    fmt(r.f_hom.gene_unq),
            'ASOC UNQ HOM F':    fmt(r.f_hom.assoc_unq),
            'DIF ASOC TOT':      dif(r.sf_hom.total,     r.f_hom.total),
            'DIF GENE UNQ':      dif(r.sf_hom.gene_unq,  r.f_hom.gene_unq),
            'DIF ASOC UNQ':      dif(r.sf_hom.assoc_unq, r.f_hom.assoc_unq),
        })
    return pd.DataFrame(data)


def build_table3(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 3: SF OG vs SF HOM."""
    data = []
    for r in rows:
        data.append({
            'Database':        r.db,
            'Relation':        r.relation,
            'GENE ID OG SF':   r.gene_id_og,
            'GENE ID HOM SF':  'ENSG',
            'ASOC ID OG':      r.assoc_id_og,
            'ASOC ID HOM':     'DOID/BTO/InChIKey',
            'ASOC OG SF':      fmt(r.sf_og.total),
            'ASOC HOM SF':     fmt(r.sf_hom.total),
            'DIF ASOC':        dif(r.sf_og.total,     r.sf_hom.total),
            'GENE UNQ OG SF':  fmt(r.sf_og.gene_unq),
            'GENE UNQ HOM SF': fmt(r.sf_hom.gene_unq),
            'DIF GENE':        dif(r.sf_og.gene_unq,  r.sf_hom.gene_unq),
            'ASOC UNQ OG SF':  fmt(r.sf_og.assoc_unq),
            'ASOC UNQ HOM SF': fmt(r.sf_hom.assoc_unq),
            'DIF ASOC UNQ':    dif(r.sf_og.assoc_unq, r.sf_hom.assoc_unq),
        })
    return pd.DataFrame(data)


def build_table4(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 4: F OG vs F HOM."""
    data = []
    for r in rows:
        data.append({
            'Database':       r.db,
            'Relation':       r.relation,
            'GENE ID OG F':   r.gene_id_og,
            'GENE ID HOM F':  'ENSG',
            'ASOC ID OG':     r.assoc_id_og,
            'ASOC ID HOM':    'DOID/BTO/InChIKey',
            'ASOC OG F':      fmt(r.f_og.total),
            'ASOC HOM F':     fmt(r.f_hom.total),
            'DIF ASOC':       dif(r.f_og.total,     r.f_hom.total),
            'GENE UNQ OG F':  fmt(r.f_og.gene_unq),
            'GENE UNQ HOM F': fmt(r.f_hom.gene_unq),
            'DIF GENE':       dif(r.f_og.gene_unq,  r.f_hom.gene_unq),
            'ASOC UNQ OG F':  fmt(r.f_og.assoc_unq),
            'ASOC UNQ HOM F': fmt(r.f_hom.assoc_unq),
            'DIF ASOC UNQ':   dif(r.f_og.assoc_unq, r.f_hom.assoc_unq),
        })
    return pd.DataFrame(data)


def build_table5(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 5: Resumen filtrado + homogenizado por database."""
    db_agg = defaultdict(lambda: {'total': 0, 'genes': set(), 'assoc': set()})
    for r in rows:
        db = r.db
        if r.f_hom.total:    db_agg[db]['total'] += r.f_hom.total
        # Approximate unique genes/assoc across relations (can't deduplicate without full data)
    data = []
    for r in rows:
        data.append({
            'Database':  r.db,
            'Relation':  r.relation,
            'ASOC TOT':  fmt(r.f_hom.total),
            'GENE UNQ':  fmt(r.f_hom.gene_unq),
            'ASOC UNQ':  fmt(r.f_hom.assoc_unq),
        })
    return pd.DataFrame(data)


def build_table6(rows: List[RelStats]) -> pd.DataFrame:
    """Tabla 6: Desglose completo por relacion, filtrado y homogenizado."""
    data = []
    for r in rows:
        data.append({
            'DATABASE':  r.db,
            'ABV':       r.relation,
            'GENE ID':   'ENSG',
            'ASOC ID':   'DOID/BTO/InChIKey',
            'ASOC TOT':  fmt(r.f_hom.total),
            'GENE UNQ':  fmt(r.f_hom.gene_unq),
            'ASOC UNQ':  fmt(r.f_hom.assoc_unq),
        })
    return pd.DataFrame(data)


# ─── Export HTML ──────────────────────────────────────────────────────────────

def df_to_html_table(df: pd.DataFrame, table_id: str) -> str:
    rows_html = ''.join(
        '<tr>' + ''.join(f'<td>{v}</td>' for v in row) + '</tr>'
        for row in df.values
    )
    headers = ''.join(f'<th>{c}</th>' for c in df.columns)
    return f'''
    <table id="{table_id}" class="display compact" style="width:100%">
      <thead><tr>{headers}</tr></thead>
      <tbody>{rows_html}</tbody>
    </table>'''


def export_html(tables: List[Tuple[str, pd.DataFrame, str]], out_path: str):
    table_ids   = [f"t{i+1}" for i in range(len(tables))]
    tab_buttons = ''.join(
        f'<button class="tablink" onclick="openTab(event,\'tab{i+1}\')">{name}</button>'
        for i, (name, _, _) in enumerate(tables)
    )
    tab_contents = ''.join(f'''
    <div id="tab{i+1}" class="tabcontent">
      <h3>{name}</h3>
      <p><em>{desc}</em></p>
      {df_to_html_table(df, f"t{i+1}")}
    </div>''' for i, (name, df, desc) in enumerate(tables))

    init_dt = ''.join(f"$('#{table_ids[i]}').DataTable({{scrollX:true,pageLength:30}});"
                      for i in range(len(tables)))

    html = f'''<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>TFM-NetActivity - Database Statistics</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; }}
  .tablink {{ background-color: #4a90d9; color: white; border: none; padding: 10px 18px;
             cursor: pointer; border-radius: 4px 4px 0 0; margin-right: 4px; }}
  .tablink:hover, .tablink.active {{ background-color: #2c5fa8; }}
  .tabcontent {{ display: none; padding: 20px; border: 1px solid #ccc; border-radius: 0 4px 4px 4px; }}
  .tabcontent.active {{ display: block; }}
  td {{ white-space: nowrap; }}
  h1 {{ color: #2c5fa8; }}
</style>
</head>
<body>
<h1>TFM-NetActivity - Database Statistics</h1>
<div style="margin-bottom:10px;">
{tab_buttons}
</div>
{tab_contents}
<script>
function openTab(evt, tabId) {{
  document.querySelectorAll('.tabcontent').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tablink').forEach(b => b.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  evt.currentTarget.classList.add('active');
}}
document.querySelector('.tablink').click();
$(document).ready(function() {{ {init_dt} }});
</script>
</body>
</html>'''
    with open(out_path, 'w') as f:
        f.write(html)
    print(f"[stats] HTML exportado -> {out_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--base-dir', default=BASE_DIR)
    p.add_argument('--out-dir',  default=None)
    return p.parse_args()


def main():
    args   = parse_args()
    base   = args.base_dir
    proc   = f"{base}/data/processed"
    out    = args.out_dir or f"{proc}/stats"
    os.makedirs(out, exist_ok=True)

    print("=" * 60)
    print(" TFM-NetActivity -- compute_stats.py")
    print("=" * 60)

    # 1. Cargar mappings
    maps = load_mappings(base)

    # 2. Compute stats per database
    all_stats: List[RelStats] = []
    all_stats += stats_opentargets(base, proc, maps)
    all_stats += stats_hpa_proteome(base, proc, maps)
    all_stats += stats_hpa_rna_cons(base, proc, maps)
    all_stats += stats_ctddisease(base, proc, maps)
    all_stats += stats_ctdchemis(base, proc, maps)
    all_stats += stats_jensentissuecurated(base, proc, maps)
    all_stats += stats_reactome(base, proc, maps)
    all_stats += stats_repohub(base, proc, maps)
    all_stats += stats_coexpressdb(base, proc, maps)
    all_stats += stats_omnipath(base, proc, maps)
    all_stats += stats_string(base, proc, maps)
    all_stats += stats_dorothea(base, proc, maps, DOROTHEA_AB, 'dorothea_AB')
    all_stats += stats_dorothea(base, proc, maps, DOROTHEA_CD, 'dorothea_CD')
    all_stats += stats_ccle_family(base, proc, maps)
    all_stats += stats_achilles(base, proc, maps)

    # 3. Build tables
    t1 = build_table1(all_stats)
    t2 = build_table2(all_stats)
    t3 = build_table3(all_stats)
    t4 = build_table4(all_stats)
    t5 = build_table5(all_stats)
    t6 = build_table6(all_stats)

    # 4. Export TSV
    for name, df in [('table1_og_sf_f', t1), ('table2_hom_sf_f', t2),
                     ('table3_comp_sf', t3), ('table4_comp_f', t4),
                     ('table5_summary', t5), ('table6_detail', t6)]:
        path = os.path.join(out, f"{name}.tsv")
        df.to_csv(path, sep='\t', index=False)
        print(f"[stats] {name}.tsv -> {path}")

    # 5. Export HTML
    tables_for_html = [
        ("T1: IDs originales SF/F",    t1, "IDs tal como estan en los archivos raw, antes y despues del filtro de calidad"),
        ("T2: IDs homogenizados SF/F",  t2, "Genes -> ENSG, Asoc -> DOID/BTO/InChIKey"),
        ("T3: Comparacion SF OG vs HOM",t3, "Perdida de informacion al homogenizar IDs (sin filtro)"),
        ("T4: Comparacion F OG vs HOM", t4, "Perdida de informacion al homogenizar IDs (filtrado)"),
        ("T5: Resumen filtrado+HOM",    t5, "Tabla final para construccion del grafo"),
        ("T6: Desglose por relacion",   t6, "Detalle por cada relacion biologica"),
    ]
    export_html(tables_for_html, os.path.join(out, "database_stats.html"))

    print("\n[stats] Completado. Archivos en:", out)


if __name__ == '__main__':
    main()
