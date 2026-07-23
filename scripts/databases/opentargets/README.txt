DATABASE:     Open Targets Platform — opentargets
VERSION:      26.06 (julio 2026)
SOURCE:       https://platform.opentargets.org
FORMAT:       Parquet (colecciones de archivos particionados)
PUBLICATIONS: https://doi.org/10.1093/nar/gkae1128

DESCARGA:
  association_overall_direct/   -> scores de asociación gen-enfermedad (direct)
  disease/                      -> metadata de enfermedades

RELATIONS:
  GEN-ass-DIS   Gene -> associates -> Disease

GENE ID:      ENSG (Ensembl Gene ID)
ASSOC ID:     EFO / MONDO ontology IDs

NOTAS TÉCNICAS:
  - Post-25.03: solo formato Parquet disponible (no JSON).
  - Se usa https:// en wget (no ftp://).
  - Requiere pandas + pyarrow para lectura.
  - association_overall_direct contiene el score armonizado final.
  - Métricas 26.03: 12,466,856 asociaciones target-enfermedad totales.
