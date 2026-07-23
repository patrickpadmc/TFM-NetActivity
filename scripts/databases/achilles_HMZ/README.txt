DATABASE:     Harmonizome (Achilles) — achilles_HMZ
SOURCE:       https://maayanlab.cloud/Harmonizome
DATASET:      Achilles Cell Line Gene Essentiality Profiles
FORMAT:       Descargado via Harmonizome API (harmonizomeapi.py)
PUBLICATION:  https://doi.org/10.1093/database/baw100

RELATIONS:    
  CLL-bfn-GEN -> Cell line -> bad fitness -> Gene
  CLL-gfn-GEN -> Cell line -> good fitness -> Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)
ASSOC ID:     Cell line name

NOTAS:
  - Mide esencialidad génica via CRISPR/RNAi knockdowns.
  - "bad fitness" = el knockdown reduce la viabilidad celular.
