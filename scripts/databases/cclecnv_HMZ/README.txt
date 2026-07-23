DATABASE:     Harmonizome (CCLE) — cclecnv_HMZ
SOURCE:       https://maayanlab.cloud/Harmonizome
DATASET:      CCLE Cell Line Gene CNV Profiles
FORMAT:       Descargado via Harmonizome API (harmonizomeapi.py)
PUBLICATION:  https://doi.org/10.1093/database/baw100

RELATIONS:
  CLL-cnd-GEN   Cell line -> copy number down -> Gene
  CLL-cnu-GEN   Cell line -> copy number up   -> Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)
ASSOC ID:     Cell line name

NOTAS:
  - up/down split se hace en script.py por el signo del CNV score.
