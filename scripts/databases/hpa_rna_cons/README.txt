DATABASE:     Human Protein Atlas — hpa_rna_cons
SOURCE:       https://www.proteinatlas.org
FILE:         rna_tissue_consensus.tsv.zip → rna_tissue_consensus.tsv
FORMAT:       TSV

RELATIONS:
  GEN-upr-TIS   Gene → upregulates → Tissue
  GEN-dwr-TIS   Gene → downregulates → Tissue

GENE ID:      ENSG (Ensembl Gene ID)
ASSOC ID:     nTPM (RNA expression consensus value)

NOTES:
  - Consenso de RNA entre múltiples fuentes del HPA.
  - up/down se determina por umbral aplicado en script.py.
