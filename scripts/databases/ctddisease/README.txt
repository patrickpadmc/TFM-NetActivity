DATABASE:     Comparative Toxigenomics Database — ctddisease
SOURCE:       https://ctdbase.org
FILE:         CTD_genes_diseases.tsv.gz
FORMAT:       TSV comprimido (gzip)

RELATIONS:
  DIS-ass-GEN   Disease → associates → Gene

GENE ID:      Entrez Gene ID (requiere mapeo a ENSG)
ASSOC ID:     Disease ID (MeSH / OMIM)

NOTES:
  - Archivo grande (~cientos de MB descomprimido). Procesado en streaming.
  - El archivo CTD_chemicals_diseases.tsv.gz lo usa ctdchemis, no este script.
