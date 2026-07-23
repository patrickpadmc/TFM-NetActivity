DATABASE:     Human Protein Atlas — hpa_proteome
SOURCE:       https://www.proteinatlas.org
FILE:         normal_ihc_data.tsv.zip -> normal_ihc_data.tsv
FORMAT:       TSV
VERSION:      25.1
PUBLICATION:  https://doi.org/10.1002/pro.4562

RELATIONS:
  GEN-pdf-TIS   Gene -> protein deficiency -> Tissue
  GEN-pab-TIS   Gene -> protein abundance  -> Tissue

GENE ID:      ENSG (Ensembl Gene ID)
ASSOC ID:     NTPM (protein expression value)

NOTES:
  - Ambas relaciones (pdf y pab) se derivan del mismo archivo fuente.
    El script.py las separa por umbral de expresión.
  - IHC: immunohistochemistry data
