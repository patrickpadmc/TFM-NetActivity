DATABASE:     Reactome — reactome
VERSION:      96
SOURCE:       https://reactome.org
FILE:         UniProt2Reactome_All_Levels.txt
FORMAT:       TSV (texto plano, sin cabecera)

RELATIONS:
  GEN-ass-PWY   Gene → associates → Pathway

GENE ID:      UniProt ID (requiere mapeo a ENSG)
ASSOC ID:     Reactome Pathway ID (R-HSA-*)

NOTAS:
  - El archivo incluye todas las especies. Filtrar por 'Homo sapiens'.
  - Columnas: UniProtID, ReactomeID, URL, description, evidenceCode, species
