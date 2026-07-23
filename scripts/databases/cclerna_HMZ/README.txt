DATABASE:     Harmonizome (CCLE) — cclerna_HMZ
SOURCE:       https://maayanlab.cloud/Harmonizome
DATASET:      CCLE Cell Line Gene Expression Profiles
FORMAT:       Descargado via Harmonizome API (harmonizomeapi.py)
DATE:         2015 APR 06
PUBLICATION:  https://doi.org/10.1093/database/baw100

RELATIONS:
  CLL-dwr-GEN   Cell line -> downregulates -> Gene
  CLL-upr-GEN   Cell line -> upregulates -> Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)
ASSOC ID:     Cell line name

NOTAS:
  - La API devuelve la URL de descarga del dataset completo.
  - El script imprime la URL real antes de descargar (útil para logs).
