# Manifiesto de reproducibilidad — graph_v4

## Identificación

- **Fecha de congelación**: 2026-07-26
- **Commit del repositorio**: `41ff89df89534d67cb382479bb53b64d36bb2500`
- **Repositorio**: `github.com/patrickpadmc/TFM-NetActivity`
- **Entorno de cómputo**: Python 3.11.3 (GCCcore-12.3.0), numpy 1.26.4, pandas 3.0.0, pyarrow 23.0.0, scipy 1.17.0, cluster SLURM UNAV (`hpclogin.unav.es`), sistema de ficheros beegfs
- **Ruta base**: `/beegfs/home/ppadmoremcc/work/TFM-NetActivity`
- **Ruta del grafo final**: `data/processed/integrated/graph_v4/graph_edges.tsv`

## Fuentes utilizadas (17 bases de datos)

| Fuente | Versión / release | Fecha del archivo fuente (verificada) | Fichero(s) principal(es) consumido(s) por script.py | SHA-256 |
|---|---|---|---|---|
| achilles_HMZ | Harmonizome, dataset original 2015-04-06 | 2026-07-25 (descarga) | gene_attribute_edges.txt.gz | 41b2791c2e05d48d5c8c2dac7b2d4467d3ddb9f99988d8549b20e32ad042e6ec |
| cclecnv_HMZ | Harmonizome, dataset original 2015-04-06 | 2026-07-25 (descarga) | gene_attribute_edges.txt.gz | 92a293f294f4a595b83b1ee74a17f36ed035f147a73993e1a1c28f607b69f532 |
| cclemut_HMZ | Harmonizome / CCLE, DepMap 26Q1 | 2026-07-25 (descarga) | OmicsSomaticMutationsMatrixDamaging.csv | 45709b1b842dd6dcd05e504cd42ca977fff760c480c60ff1184c1a5935b595b0 |
| | | | OmicsSomaticMutationsMatrixHotspot.csv | a86cd8c92b86d5507e63103f06a507913897c85262e0a0d7e44fa22768dd4dc9 |
| cclerna_HMZ | Harmonizome / CCLE, DepMap 26Q1 | 2026-07-25 (descarga) | OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | 0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a |
| coexpressdb | COXPRESdb v22-05 (Hsa-r.c6-0) | 2026-05-22 (descarga) | Hsa-r.v22-05.zip | 2533d8cf9dc30bbb1d60a188b3a4a2b52872d11eaf898a1a568be1c743e9a5a5 |
| ctdchemgen | CTD, CTD_chem_gene_ixns | 2026-05-28 (descarga) | CTD_chem_gene_ixns.tsv.gz | cfb8091e6d36396e72ba3aeba32f66aec777af5a47c73ff0ada9f4a82e65b9f1 |
| ctdchemis | CTD, CTD_chemicals_diseases | 2026-06-29 (descarga) | CTD_chemicals_diseases.tsv.gz | ff95b6c15c774b1168904b61ae07e7d4e9b0d7ede8df049d44fe740a58494c27 |
| dorothea_AB | DoRothEA via OmniPath, niveles A/B | 2026-06-25 (símlink a omnipath) | omnipath_webservice_interactions__latest.tsv.gz | e1f2777f5ac7c10390402c442b54621b6fc3de761b531ab91b46330516875f3a |
| dorothea_CD | DoRothEA via OmniPath, niveles C/D | 2026-06-25 (símlink a omnipath) | omnipath_webservice_interactions__latest.tsv.gz | e1f2777f5ac7c10390402c442b54621b6fc3de761b531ab91b46330516875f3a |
| hpa_proteome | Human Protein Atlas v25.1 | verificado 2026-06-02, contenido idéntico al servido actualmente | normal_ihc_data.tsv | 55415dd27fe812242d43add4ee7b24ec6b0514a314c38c147b2813d23cb4b801 |
| hpa_rna_cons | Human Protein Atlas v25.1 | verificado 2025-11-06, contenido idéntico al servido actualmente (comprobado vs. Content-Length/Last-Modified remoto) | rna_tissue_consensus.tsv | cdedaeaf3cdfc89e22b3891ea24ae2afabc0afd26d8883076121a363608450b6 |
| jensentissuecurated | Jensen Lab COMPARTMENTS, canal knowledge | 2026-07-10 (descarga) | human_tissue_knowledge_full.tsv | fd92383a27cd1df858bfd8c78f0e354fa5cea5ceda5b03c81b4633b937af93c5 |
| omnipath | OmniPath (archive.omnipathdb.org) | interactions: 2026-06-25; enz_sub: verificado 2025-08-13, contenido idéntico al servido actualmente | omnipath_webservice_interactions__latest.tsv.gz | e1f2777f5ac7c10390402c442b54621b6fc3de761b531ab91b46330516875f3a |
| | | | omnipath_webservice_enz_sub__latest.tsv.gz | 10c889259078136e5572a75783343000998e33040f7ed3037aed4d5a9ecab32a |
| opentargets | Open Targets Platform 26.06 | 2026-06-22 (descarga, 39 ficheros part-*.parquet) | association_overall_indirect/part-*.parquet (39 ficheros) | (39 hashes individuales, no listados por brevedad; ver `sha256sum data/raw/databases/opentargets/part-*.parquet` en el cluster) |
| reactome | Reactome v97 | 2026-06-24 (descarga) | UniProt2Reactome_All_Levels.txt | a9c41c9b997637286a06b9acf69557e8c2c6fc703599e3120dd6a85bfe3736f7 |
| repohub | The Drug Repurposing Hub (dataset estático, 2020-03-24) | 2026-07-25 (subida manual, sin versión más nueva disponible) | repo-drug-annotation-20200324.txt | 5f8284538e73c19a316d1cf45ea10300de6e15c555d75ced7c1b7387457fc523 |
| string | STRING v12.0 | verificado 2023-05-16, contenido idéntico al servido actualmente (sin release más nueva) | 9606.protein.links.full.v12.0.txt.gz | 07f9fa42ae5006ccb0b4694c17b97dd4d346c141bf769006cb6ffc4c8e04d016 |

*Nota sobre "fecha del archivo fuente": para string, hpa_rna_cons y omnipath/enz_sub, la fecha no corresponde al momento de la descarga sino a la última actualización de la fuente original, se verificó explícitamente (vía `curl -I`, comparando `Content-Length` y `Last-Modified`) que el contenido local es idéntico, byte a byte, al que sirve la fuente hoy.*

## Nodos por tipo

| Tipo | N |
|---|---|
| GEN | 26,330 |
| DIS | 23,697 |
| CPD | 19,629 |
| PWY | 2,820 |
| CLL | 1,953 |
| TIS | 699 |
| **Total** | **75,128** |

## Aristas por relación canónica

Ver tabla completa (38 filas) en `seccion1_auditoria_graph_v4.md`, Sección 1. Total: **110,690,581** aristas, **34** tipos de relación canónicos, **17** bases de datos.

## Reglas de filtrado y deduplicación

- Cada base de datos aplica sus propios filtros de calidad de forma independiente, dentro de su `script.py`, alineados a la metodología de Bioteque (IRB Barcelona) y no se aplica ningún filtro global sobre el grafo ya integrado.
- Deduplicación por par (node1, node2) dentro de cada combinación fuente+relación, verificada empíricamente tras la integración: 0 duplicados exactos de fila completa, 0 duplicados de par dentro de la misma fuente+relación, sobre las 110,690,581 aristas totales.
- Identificadores homogenizados por tipo de nodo: GEN -> ENSG (Ensembl Gene ID), CPD -> InChIKey (con fallback a ID crudo de CTD cuando no hay mapeo disponible en Bioteque (ver limitación conocida más abajo)), CLL -> CVCL/RRID (DepMap), TIS -> BTO, DIS -> ontología nativa de cada fuente (MONDO/EFO/MeSH).
- `coexpressdb` (GEN-cex-GEN): pares simetrizados mediante `tuple(sorted(...))`; en colisiones de pares duplicados (evidencia en ambos sentidos), el valor de `mutual_rank` se combina como `sqrt(a * b)` (media geométrica), con los valores negativos de origen (artefacto de piso del score log-transformado de COXPRESdb) fijados a 0.0 antes de combinar.

## Limitaciones conocidas (no bloqueantes, documentadas)

- Cobertura de compuestos: no todos los compuestos de CTD/repohub tienen mapeo a InChIKey disponible en Bioteque; los que no lo tienen retienen su identificador crudo de CTD (formato MeSH, prefijo C o D indistintamente. El prefijo no distingue compuesto de enfermedad en MeSH, solo distingue encabezado principal de registro suplementario).
- 28 nodos (14 pares) quedan fuera del componente conexo principal (99.96% del grafo), por escasez genuina de datos (evidencia única, sin corroboración cruzada entre fuentes), no por error de homogenización. Documentados para excluir de los splits de evaluación de predicción de enlaces.

## Semilla global

`seed=42` : convención establecida del proyecto para todo paso con aleatoriedad (muestreo, splits de train/test, inicialización de modelos). La construcción de graph_v4 en sí (concatenación de las 38 tablas ya procesadas) es determinista y no usa aleatoriedad.

## Scripts de construcción y verificación

- `scripts/integrated/graph_v4/build_graph_table_v4.py` — integración de las 38 tablas fuente en `graph_edges.tsv`.
- `scripts/analysis/connected_components_v4.py` — análisis de componentes conexas.
- `scripts/analysis/audit_graph_v4.py` — auditoría de duplicados, solapamiento semántico, direccionalidad empírica y nodos únicos por relación.
