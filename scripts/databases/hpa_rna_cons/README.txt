DATABASE:     Human Protein Atlas — hpa_rna_cons
SOURCE:       https://www.proteinatlas.org
FILE:         rna_tissue_consensus.tsv.zip -> rna_tissue_consensus.tsv
FORMAT:       TSV
VERSION:      25.1
PUBLICATION:  https://doi.org/10.1002/pro.4562

RELATIONS:
  GEN-upr-TIS   Gene -> upregulates -> Tissue
  GEN-dwr-TIS   Gene -> downregulates -> Tissue

GENE ID:      ENSG (Ensembl Gene ID)
ASSOC ID:     nTPM (RNA expression consensus value)

NOTES:
  - Consenso de RNA entre múltiples fuentes del HPA.
  - up/down se determina por umbral aplicado en script.py.
  - Primero se escala cada vector de expresion de genes a traves de los tejidos y luego se seleccionan los mejores 250 upd y down-regulated genes para cada tejido. La expresion del Z-score entre -0.5 y 0.5 fueron saltados.
