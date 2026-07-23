DATABASE:     STRING — string
VERSION:      v12.0
SOURCE:       https://string-db.org
ARCHIVOS:
  9606.protein.links.full.v12.0.txt.gz   -> interacciones PPI con todos los scores
  human.uniprot_2_string.2018.tsv.gz     -> mapeo UniProt -> STRING protein IDs
DATE:         2026 JUL
PUBLICATION:  https://doi.org/10.1093/nar/gkac1000

RELATIONS:
  GEN-ppi-GEN   Gene -> protein-protein interaction -> Gene

GENE ID:      ENSG (tras mapeo STRING ID -> UniProt -> ENSG)
ASSOC ID:     combined_score (0–1000)

NOTAS:
  - El archivo full incluye scores por canal: cooccurence, coexpression,
    experimental, database, textmining, etc.
  - Filtrar por taxon 9606 (Homo sapiens) si el archivo incluye otras especies.
  - El mapeo de 2018 es el más reciente disponible en string-db.org para UniProt.
