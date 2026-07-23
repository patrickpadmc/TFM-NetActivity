#!/usr/bin/env python3
"""
ctdchemgen - script.py
Exploracion de CTD chemical-gene interactions (CTD_chem_gene_ixns.tsv.gz).

Esta NO es la version final de procesado. Su unico objetivo es responder
a la pregunta: que tipos de interaccion (campo InteractionActions) hay
en el fichero y con que frecuencia, para decidir cuales son las
puntuales a integrar en la KG (relaciones CPD-???-GEN).

Columnas del fichero (tras las lineas de comentario que empiezan por '#'):
  ChemicalName | ChemicalID | CasRN | GeneSymbol | GeneID | GeneForms |
  Organism | OrganismID | Interaction | InteractionActions | PubMedIDs

InteractionActions es una lista separada por '|' de pares "degree^type"
(ej. "increases^expression|affects^binding"). degree es habitualmente
increases / decreases / affects.

Este script NO usa pandas/numpy: lee el gzip en streaming linea a linea,
por lo que puede ejecutarse en el nodo de login sin problema (no
requiere srun). Si una version posterior anade pandas/numpy para el
filtrado final, esa parte si debera lanzarse via srun.

Output (en --out-dir):
  summary.tsv                    -> conteos generales
  organism_counts.tsv            -> filas por organismo
  degree_counts.tsv              -> filas por degree (increases/decreases/affects)
  base_type_counts.tsv           -> filas por tipo base (expression, binding, ...)
  interactionactions_counts.tsv  -> filas por combinacion completa degree^type

Uso:
  python script.py --raw-dir <ruta_raw> --out-dir <ruta_out>
                    [--filename CTD_chem_gene_ixns.tsv.gz]
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
IDX_ORGANISM = COL_NAMES.index("Organism")
IDX_INTERACTION_ACTIONS = COL_NAMES.index("InteractionActions")
IDX_CHEMICAL_ID = COL_NAMES.index("ChemicalID")
IDX_GENE_ID = COL_NAMES.index("GeneID")


def parse_args():
    p = argparse.ArgumentParser(
        description="Explora los tipos de interaccion presentes en CTD_chem_gene_ixns.tsv.gz"
    )
    p.add_argument("--raw-dir", required=True, help="Carpeta donde esta el tsv.gz descargado")
    p.add_argument("--out-dir", required=True, help="Carpeta donde se escriben las tablas de stats")
    p.add_argument(
        "--filename",
        default="CTD_chem_gene_ixns.tsv.gz",
        help="Nombre del fichero descargado dentro de --raw-dir",
    )
    return p.parse_args()


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
            "ERROR: %s no es un gzip valido (no tiene la cabecera gzip). Es "
            "posible que la descarga haya devuelto una pagina de verificacion "
            "(captcha) en vez del fichero real. Borra el fichero y descargalo "
            "manualmente desde el navegador, luego subelo con scp.\n" % path
        )
        sys.exit(1)
    return gzip.open(path, "rt", encoding="utf-8", errors="replace")


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    in_path = os.path.join(args.raw_dir, args.filename)

    total_lines = 0
    data_rows = 0
    bad_ncols = 0
    empty_actions = 0

    organism_counts = Counter()
    degree_counts = Counter()
    base_type_counts = Counter()
    full_action_counts = Counter()

    chemicals = set()
    genes = set()

    with open_input(in_path) as f:
        for line in f:
            total_lines += 1
            if line.startswith("#") or not line.strip():
                continue
            data_rows += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) != EXPECTED_NCOLS:
                bad_ncols += 1
                continue

            organism_counts[cols[IDX_ORGANISM]] += 1
            chemicals.add(cols[IDX_CHEMICAL_ID])
            genes.add(cols[IDX_GENE_ID])

            actions = cols[IDX_INTERACTION_ACTIONS].strip()
            if not actions:
                empty_actions += 1
                continue

            for piece in actions.split("|"):
                piece = piece.strip()
                if not piece:
                    continue
                full_action_counts[piece] += 1
                if "^" in piece:
                    degree, base_type = piece.split("^", 1)
                else:
                    degree, base_type = "unspecified", piece
                degree_counts[degree] += 1
                base_type_counts[base_type] += 1

            if data_rows % 500000 == 0:
                sys.stderr.write("...%d filas procesadas\n" % data_rows)

    def write_counter(counter, path, header):
        with open(path, "w") as o:
            o.write("%s\tcount\n" % header)
            for key, count in counter.most_common():
                o.write("%s\t%d\n" % (key, count))

    write_counter(organism_counts, os.path.join(args.out_dir, "organism_counts.tsv"), "organism")
    write_counter(degree_counts, os.path.join(args.out_dir, "degree_counts.tsv"), "degree")
    write_counter(base_type_counts, os.path.join(args.out_dir, "base_type_counts.tsv"), "base_type")
    write_counter(
        full_action_counts,
        os.path.join(args.out_dir, "interactionactions_counts.tsv"),
        "interaction_action",
    )

    with open(os.path.join(args.out_dir, "summary.tsv"), "w") as o:
        o.write("metric\tvalue\n")
        o.write("total_lines_in_file\t%d\n" % total_lines)
        o.write("data_rows\t%d\n" % data_rows)
        o.write("rows_bad_ncols\t%d\n" % bad_ncols)
        o.write("rows_empty_interactionactions\t%d\n" % empty_actions)
        o.write("unique_chemicals\t%d\n" % len(chemicals))
        o.write("unique_genes\t%d\n" % len(genes))
        o.write("unique_organisms\t%d\n" % len(organism_counts))
        o.write("unique_base_types\t%d\n" % len(base_type_counts))
        o.write("unique_full_actions\t%d\n" % len(full_action_counts))

    sys.stderr.write("Done. Stats escritas en %s\n" % args.out_dir)


if __name__ == "__main__":
    main()
