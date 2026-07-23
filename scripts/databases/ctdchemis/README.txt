DATABASE:     Comparative Toxigenomics Database — ctdchemis
SOURCE:       https://ctdbase.org
FILE:         CTD_chemicals_diseases.tsv.gz
FORMAT:       TSV comprimido (gzip)

RELATIONS:
  CPD-cau-DIS   Compound → causes  → Disease
  CPD-trt-DIS   Compound → treats  → Disease

GENE ID:      N/A (relación compuesto–enfermedad)
ASSOC ID:     Disease ID (MeSH), Chemical ID (MeSH)

NOTES:
  - La separación causes/treats se hace por el campo DirectEvidence en script.py.
  - Archivo potencialmente muy grande; procesado en streaming.
