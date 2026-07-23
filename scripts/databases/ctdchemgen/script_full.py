#!/usr/bin/env python3
"""
ctdchemgen - script.py
Procesado de CTD chemical-gene interactions (CTD_chem_gene_ixns.tsv.gz).

Logica:
  1. Filtrar a Homo sapiens (OrganismID 9606), salvo --all-organisms.
  2. Mapear:
       ChemicalID (CTD, formato MESH Cxxxxx/Dxxxxx) -> InChIKey
         via bioteque metadata/mappings/CPD/ctd.tsv
         fallback: se mantiene el ID de CTD si no hay InChIKey (igual que
         en ctdchemis).
       GeneID (Entrez) -> ENSG
         via data/processed/entrez2ensg.tsv
         sin fallback: si no hay ENSG, la fila se descarta (el universo
         de genes del proyecto es ENSG humano).
  3. Parsear InteractionActions ("degree^type", lista separada por '|')
     y mapear cada (degree, type) a una relacion CPD-???-GEN segun la
     tabla GROUPS / RELATIONS de abajo.

Se descartan estos tipos (no aportan una relacion binaria limpia o
tienen volumen marginal/ambiguo): cotreatment, response to substance,
splicing, stability, folding, mutagenesis.

Columnas del fichero de entrada (tras las lineas de comentario '#'):
  ChemicalName | ChemicalID | CasRN | GeneSymbol | GeneID | GeneForms |
  Organism | OrganismID | Interaction | InteractionActions | PubMedIDs

Output (en --out-dir), un fichero por relacion con columnas n1 (InChIKey
o ID de CTD si no mapeo) n2 (ENSG):
  CPD-upe-GEN.tsv  increases^expression
  CPD-dwe-GEN.tsv  decreases^expression
  CPD-afe-GEN.tsv  affects^expression
  CPD-upa-GEN.tsv  increases^activity
  CPD-dwa-GEN.tsv  decreases^activity
  CPD-afa-GEN.tsv  affects^activity
  CPD-upb-GEN.tsv  increases^abundance
  CPD-dwb-GEN.tsv  decreases^abundance
  CPD-afb-GEN.tsv  affects^abundance
  CPD-bnd-GEN.tsv  binding (cualquier degree)
  CPD-upm-GEN.tsv  increases^{reaction, metabolic processing, o cualquier
                   modificacion quimica especifica: methylation,
                   phosphorylation, oxidation, acetylation, ubiquitination,
                   sumoylation, glycosylation*, glucuronidation,
                   palmitoylation, sulfation, nitrosation, lipidation,
                   prenylation, farnesylation, geranoylation,
                   carboxylation, carbamoylation, amination, alkylation,
                   acylation, ethylation, ribosylation, ADP-ribosylation,
                   glutathionylation, reduction, hydrolysis,
                   hydroxylation, cleavage, degradation, chemical
                   synthesis, glycation, myristoylation, polymerization}
  CPD-dwm-GEN.tsv  decreases^{idem}
  CPD-afm-GEN.tsv  affects^{idem}
  CPD-loc-GEN.tsv  {localization, transport, secretion, uptake, export,
                   import} (cualquier degree, relacion sin direccion)
  stats.tsv        filas escritas por relacion + descartadas + sin mapear

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_out>
                    [--filename CTD_chem_gene_ixns.tsv.gz]
                    [--ctd-mapping <ruta_a_CPD/ctd.tsv>]
                    [--entrez-mapping <ruta_a_entrez2ensg.tsv>]
                    [--all-organisms]
"""
import os
import sys
import gzip
import argparse
from collections import Counter

EXPECTED_NCOLS = 11
COL_NAMES = [
    "ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
    "GeneForms", "Organism", "OrganismID", "Interaction",
    "InteractionActions", "PubMedIDs",
]
IDX_CHEMICAL_ID = COL_NAMES.index("ChemicalID")
IDX_GENE_ID = COL_NAMES.index("GeneID")
IDX_ORGANISM_ID = COL_NAMES.index("OrganismID")
IDX_INTERACTION_ACTIONS = COL_NAMES.index("InteractionActions")

HUMAN_TAXID = "9606"

# base_type (tal cual aparece en InteractionActions) -> grupo
BASE_TYPE_TO_GROUP = {
    "expression": "expression",
    "activity": "activity",
    "abundance": "abundance",
    "binding": "binding",
    "reaction": "biotransformation",
    "metabolic processing": "biotransformation",
    "methylation": "biotransformation",
    "phosphorylation": "biotransformation",
    "oxidation": "biotransformation",
    "acetylation": "biotransformation",
    "ubiquitination": "biotransformation",
    "sumoylation": "biotransformation",
    "glycosylation": "biotransformation",
    "O-linked glycosylation": "biotransformation",
    "N-linked glycosylation": "biotransformation",
    "glucuronidation": "biotransformation",
    "palmitoylation": "biotransformation",
    "sulfation": "biotransformation",
    "nitrosation": "biotransformation",
    "lipidation": "biotransformation",
    "prenylation": "biotransformation",
    "farnesylation": "biotransformation",
    "geranoylation": "biotransformation",
    "carboxylation": "biotransformation",
    "carbamoylation": "biotransformation",
    "amination": "biotransformation",
    "alkylation": "biotransformation",
    "acylation": "biotransformation",
    "ethylation": "biotransformation",
    "ribosylation": "biotransformation",
    "ADP-ribosylation": "biotransformation",
    "glutathionylation": "biotransformation",
    "reduction": "biotransformation",
    "hydrolysis": "biotransformation",
    "hydroxylation": "biotransformation",
    "cleavage": "biotransformation",
    "degradation": "biotransformation",
    "chemical synthesis": "biotransformation",
    "glycation": "biotransformation",
    "myristoylation": "biotransformation",
    "polymerization": "biotransformation",
    "localization": "trafficking",
    "transport": "trafficking",
    "secretion": "trafficking",
    "uptake": "trafficking",
    "export": "trafficking",
    "import": "trafficking",
}
# Descartados a proposito (no quedan en BASE_TYPE_TO_GROUP, se ignoran):
#   cotreatment, response to substance, splicing, stability, folding,
#   mutagenesis

DIRECTED_GROUPS = {
    "expression": {"increases": "CPD-upe-GEN", "decreases": "CPD-dwe-GEN", "affects": "CPD-afe-GEN"},
    "activity": {"increases": "CPD-upa-GEN", "decreases": "CPD-dwa-GEN", "affects": "CPD-afa-GEN"},
    "abundance": {"increases": "CPD-upb-GEN", "decreases": "CPD-dwb-GEN", "affects": "CPD-afb-GEN"},
    "biotransformation": {"increases": "CPD-upm-GEN", "decreases": "CPD-dwm-GEN", "affects": "CPD-afm-GEN"},
}
UNDIRECTED_GROUPS = {
    "binding": "CPD-bnd-GEN",
    "trafficking": "CPD-loc-GEN",
}

ALL_RELATIONS = sorted(
    set(v for g in DIRECTED_GROUPS.values() for v in g.values())
    | set(UNDIRECTED_GROUPS.values())
)


def parse_args():
    p = argparse.ArgumentParser(description="Procesa CTD_chem_gene_ixns.tsv.gz a relaciones CPD-???-GEN")
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--filename", default="CTD_chem_gene_ixns.tsv.gz")
    p.add_argument(
        "--ctd-mapping",
        default=os.path.expanduser(
            "~/work/external/bioteque/metadata/mappings/CPD/ctd.tsv"
        ),
        help="Ruta a CPD/ctd.tsv (CTD ChemicalID -> InChIKey) de Bioteque",
    )
    p.add_argument(
        "--entrez-mapping",
        default=os.path.expanduser(
            "~/work/TFM-NetActivity/data/processed/entrez2ensg.tsv"
        ),
        help="Ruta a entrez2ensg.tsv",
    )
    p.add_argument("--all-organisms", action="store_true", help="No filtrar por OrganismID 9606")
    return p.parse_args()


def load_ctd2ikey(path):
    d = {}
    if not os.path.isfile(path):
        sys.stderr.write("AVISO: no se encontro %s, no habra mapeo a InChIKey (fallback a ID CTD)\n" % path)
        return d
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3 or not cols[2]:
                continue
            d[cols[0]] = cols[2]
    return d


def load_entrez2ensg(path):
    d = {}
    if not os.path.isfile(path):
        sys.stderr.write("ERROR: no se encontro %s (necesario para mapear genes a ENSG)\n" % path)
        sys.exit(1)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2 or not cols[1]:
                continue
            d[cols[0]] = cols[1]
    return d


def open_input(path):
    if not os.path.isfile(path):
        sys.stderr.write("ERROR: no existe %s\n" % path)
        sys.exit(1)
    if os.path.getsize(path) == 0:
        sys.stderr.write("ERROR: %s esta vacio (descarga incompleta?)\n" % path)
        sys.exit(1)
    with open(path, "rb") as raw:
        magic = raw.read(2)
    if magic != b"\x1f\x8b":
        sys.stderr.write(
            "ERROR: %s no es un gzip valido (posible pagina de captcha en vez "
            "del fichero real). Borra el fichero y descargalo manualmente.\n" % path
        )
        sys.exit(1)
    return gzip.open(path, "rt", encoding="utf-8", errors="replace")


def resolve_relation(base_type, degree):
    group = BASE_TYPE_TO_GROUP.get(base_type)
    if group is None:
        return None
    if group in UNDIRECTED_GROUPS:
        return UNDIRECTED_GROUPS[group]
    by_degree = DIRECTED_GROUPS[group]
    return by_degree.get(degree, by_degree["affects"])


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    in_path = os.path.join(args.raw_dir, args.filename)

    ctd2ikey = load_ctd2ikey(args.ctd_mapping)
    entrez2ensg = load_entrez2ensg(args.entrez_mapping)

    pairs = {rel: set() for rel in ALL_RELATIONS}

    n_total = 0
    n_skipped_organism = 0
    n_skipped_bad_ncols = 0
    n_skipped_no_gene_map = 0
    n_skipped_no_action = 0
    n_actions_discarded_type = 0
    n_chemicals_unmapped = set()

    with open_input(in_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            n_total += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) != EXPECTED_NCOLS:
                n_skipped_bad_ncols += 1
                continue

            if not args.all_organisms and cols[IDX_ORGANISM_ID] != HUMAN_TAXID:
                n_skipped_organism += 1
                continue

            gene_id = cols[IDX_GENE_ID]
            ensg = entrez2ensg.get(gene_id)
            if not ensg:
                n_skipped_no_gene_map += 1
                continue

            chem_id = cols[IDX_CHEMICAL_ID]
            if chem_id in ctd2ikey:
                n1 = ctd2ikey[chem_id]
            else:
                n1 = chem_id
                n_chemicals_unmapped.add(chem_id)

            actions = cols[IDX_INTERACTION_ACTIONS].strip()
            if not actions:
                n_skipped_no_action += 1
                continue

            for piece in actions.split("|"):
                piece = piece.strip()
                if not piece or "^" not in piece:
                    continue
                degree, base_type = piece.split("^", 1)
                relation = resolve_relation(base_type, degree)
                if relation is None:
                    n_actions_discarded_type += 1
                    continue
                pairs[relation].add((n1, ensg))

            if n_total % 500000 == 0:
                sys.stderr.write("...%d filas leidas\n" % n_total)

    for rel in ALL_RELATIONS:
        out_path = os.path.join(args.out_dir, "%s.tsv" % rel)
        with open(out_path, "w") as o:
            o.write("n1\tn2\n")
            for n1, n2 in sorted(pairs[rel]):
                o.write("%s\t%s\n" % (n1, n2))

    with open(os.path.join(args.out_dir, "stats.tsv"), "w") as o:
        o.write("metric\tvalue\n")
        o.write("rows_read\t%d\n" % n_total)
        o.write("rows_skipped_bad_ncols\t%d\n" % n_skipped_bad_ncols)
        o.write("rows_skipped_non_human\t%d\n" % n_skipped_organism)
        o.write("rows_skipped_no_gene_mapping\t%d\n" % n_skipped_no_gene_map)
        o.write("rows_skipped_no_interactionactions\t%d\n" % n_skipped_no_action)
        o.write("actions_discarded_unmapped_type\t%d\n" % n_actions_discarded_type)
        o.write("chemicals_without_inchikey_kept_as_ctdid\t%d\n" % len(n_chemicals_unmapped))
        for rel in ALL_RELATIONS:
            o.write("pairs_%s\t%d\n" % (rel, len(pairs[rel])))

    sys.stderr.write("Done. Output en %s\n" % args.out_dir)


if __name__ == "__main__":
    main()
