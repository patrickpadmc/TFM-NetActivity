DATABASE:     DoRothEA via OmniPath — dorothea_AB
SOURCE:       https://omnipathdb.org / https://archive.omnipathdb.org
FILE:         omnipath_webservice_interactions__latest.tsv.gz
              (symlink al archivo descargado en omnipath/)
FORMAT:       TSV comprimido (gzip)
DATE:         JUN 25 2026
PUBLICATION:  https://doi.org/10.1101/gr.240663.118

RELATIONS:
  GEN-upr-GEN   Gene -> upregulates   -> Gene  (TF -> target)
  GEN-reg-GEN   Gene -> regulates     -> Gene
  GEN-dwr-GEN   Gene -> downregulates -> Gene

GENE ID:      Gene symbol (requiere mapeo a ENSG)
CONFIDENCE:   Niveles A y B (alta confianza)

NOTAS:
  - DoRothEA está integrado dentro de OmniPath. El filtrado por nivel A/B
    se hace en script.py usando la columna 'dorothea_level'.
  - El archivo fuente es compartido con omnipath/ (symlink, no duplicado).
  - A: ChIP-seq experimental + literatura curada
  - B: ChIP-seq o literatura curada
