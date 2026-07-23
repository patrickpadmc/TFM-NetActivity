DATABASE:     DoRothEA via OmniPath — dorothea_CD
SOURCE:       https://omnipathdb.org / https://archive.omnipathdb.org
FILE:         omnipath_webservice_interactions__latest.tsv.gz
              (symlink al archivo descargado en omnipath/)
FORMAT:       TSV comprimido (gzip)

RELATIONS:
  GEN-upr-GEN   Gene → upregulates  → Gene  (TF → target)
  GEN-reg-GEN   Gene → regulates    → Gene
  GEN-dwr-GEN   Gene → downregulates→ Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)
CONFIDENCE:   Niveles C y D (confianza media-baja)

NOTAS:
  - Igual que dorothea_AB pero filtrando por dorothea_level in ['C', 'D'].
  - C: solo predicciones computacionales curadas
  - D: solo predicciones computacionales
