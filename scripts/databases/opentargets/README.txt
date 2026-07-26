DATABASE:     Open Targets Platform — opentargets
VERSION:      26.06 (julio 2026)
SOURCE:       https://platform.opentargets.org
FORMAT:       Parquet (colecciones de archivos particionados)
PUBLICATIONS: https://doi.org/10.1093/nar/gkae1128

DESCARGA:
  association_overall_indirect/ -> scores de asociación gen-enfermedad (indirect, incluye inferencias via EFO/MONDO)
  disease/                      -> metadata de enfermedades

RELATIONS:
  GEN-ass-DIS   Gene -> associates -> Disease

GENE ID:      ENSG (Ensembl Gene ID)
ASSOC ID:     EFO / MONDO ontology IDs

NOTAS TÉCNICAS:
  - Post-25.03: solo formato Parquet disponible (no JSON).
  - Se usa ftp:// en wget (https:// se cuelga tras el TLS handshake en este cluster, nunca responde).
  - Requiere pandas + pyarrow para lectura.
  - association_overall_indirect contiene el score armonizado final (incluye asociaciones inferidas por propagacion ontologica via EFO/MONDO).
  - Métricas 26.03: 12,466,856 asociaciones target-enfermedad totales.
