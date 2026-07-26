# Sección 1 — Auditoría y congelación de graph_v4

**Fecha de auditoría**: 2026-07-26
**Commit del repositorio**: `41ff89df89534d67cb382479bb53b64d36bb2500`
**Repositorio**: `github.com/patrickpadmc/TFM-NetActivity`

---

## 1. Tabla de auditoría de relaciones

38 combinaciones fuente+relación, verificadas directamente contra el archivo final `graph_edges.tsv` (110,690,581 filas) y contra los scripts de procesamiento (`script.py`) de cada base de datos. Ningún número de esta tabla es estimado: todos provienen de conteos ejecutados sobre los datos reales.

| Fuente | Relación original (fuente) | Relación canónica | Nodo origen | Nodo destino | Aristas | Nodos únicos | Dirección | Estado |
|---|---|---|---|---|---|---|---|---|
| opentargets | association_overall_indirect (score armonizado gen-enfermedad, incl. inferencias EFO/MONDO) | GEN-ass-DIS | GEN | DIS | 283,977 | 31,624 | Dirigida (bipartita) | Integrado |
| hpa_proteome | Abundancia de proteína por IHC (normal_ihc_data.tsv) | GEN-pab-TIS | GEN | TIS | 77,568 | 8,632 | Dirigida (bipartita) | Integrado |
| hpa_proteome | Deficiencia de proteína por IHC | GEN-pdf-TIS | GEN | TIS | 103,702 | 9,654 | Dirigida (bipartita) | Integrado |
| hpa_rna_cons | Top 250 genes sobreexpresados por tejido, consenso RNA (z-score) | GEN-upr-TIS | GEN | TIS | 12,750 | 5,389 | Dirigida (bipartita) | Integrado |
| hpa_rna_cons | Top 250 genes subexpresados por tejido, consenso RNA | GEN-dwr-TIS | GEN | TIS | 12,750 | 4,632 | Dirigida (bipartita) | Integrado |
| jensentissuecurated | Canal "knowledge" (curado por literatura/bases de datos) | GEN-ass-TIS | GEN | TIS | 315,368 | 17,277 | Dirigida (bipartita) | Integrado |
| reactome | UniProt2Reactome_All_Levels.txt, filtrado Homo sapiens | GEN-ass-PWY | GEN | PWY | 152,359 | 15,633 | Dirigida (bipartita) | Integrado |
| repohub | Drug Repurposing Hub, target de fármaco | CPD-int-GEN | CPD | GEN | 10,279 | 6,492 | Dirigida (bipartita) | Integrado |
| coexpressdb | Mutual rank de coexpresión (COXPRESdb v22-05) | GEN-cex-GEN | GEN | GEN | 105,907,896 | 14,557 | No dirigida (simetrizada por construcción, par ordenado) | Integrado |
| omnipath | Interacción proteína-proteína | GEN-ppi-GEN | GEN | GEN | 84,608 | 13,557 | No dirigida semánticamente; almacenada en un solo sentido por par (0% invertidos) | Integrado — mantenido separado de string (ver Decisión 1) |
| omnipath | Fosforilación (interactions) | GEN-pho-GEN | GEN | GEN | 27,639 | 5,427 | Dirigida (3.83% de pares bidireccionales, biológicamente plausible) | Integrado |
| omnipath | Desfosforilación (enzyme-substrate) | GEN-dph-GEN | GEN | GEN | 749 | 395 | Dirigida (0.27% bidireccional) | Integrado |
| string | Interacción proteína-proteína, combined_score | GEN-ppi-GEN | GEN | GEN | 229,073 | 15,835 | No dirigida semánticamente; almacenada en un solo sentido por par (0% invertidos) | Integrado — mantenido separado de omnipath (ver Decisión 1) |
| dorothea_AB | Regulón TF->gen, confianza A/B (interactions, campo dorothea_level) | GEN-reg-GEN | GEN | GEN | 16,355 | 5,741 | Dirigida (0.83% bidireccional) | Integrado — mantenido separado de dorothea_CD (ver Decisión 2) |
| dorothea_AB | Activación TF->gen, confianza A/B | GEN-upr-GEN | GEN | GEN | 6,528 | 3,240 | Dirigida (0% bidireccional) | Integrado — mantenido separado de dorothea_CD |
| dorothea_AB | Represión TF->gen, confianza A/B | GEN-dwr-GEN | GEN | GEN | 1,097 | 888 | Dirigida (0% bidireccional) | Integrado — mantenido separado de dorothea_CD |
| dorothea_CD | Regulón TF->gen, confianza C/D | GEN-reg-GEN | GEN | GEN | 286,551 | 20,031 | Dirigida (0.43% bidireccional) | Integrado — mantenido separado de dorothea_AB |
| dorothea_CD | Activación TF->gen, confianza C/D | GEN-upr-GEN | GEN | GEN | 12,640 | 5,230 | Dirigida (0.47% bidireccional) | Integrado — mantenido separado de dorothea_AB |
| dorothea_CD | Represión TF->gen, confianza C/D | GEN-dwr-GEN | GEN | GEN | 3,505 | 2,064 | Dirigida (0.11% bidireccional) | Integrado — mantenido separado de dorothea_AB |
| cclerna_HMZ | CCLE, gen sobreexpresado en línea celular | CLL-upr-GEN | CLL | GEN | 401,250 | 6,905 | Dirigida (bipartita) | Integrado |
| cclerna_HMZ | CCLE, gen subexpresado en línea celular | CLL-dwr-GEN | CLL | GEN | 401,250 | 11,199 | Dirigida (bipartita) | Integrado |
| cclemut_HMZ | CCLE, mutación Hotspot/Damaging | CLL-mut-GEN | CLL | GEN | 107,951 | 19,985 | Dirigida (bipartita) | Integrado |
| cclecnv_HMZ | CCLE, ganancia de número de copias | CLL-cnu-GEN | CLL | GEN | 508,053 | 21,874 | Dirigida (bipartita) | Integrado |
| cclecnv_HMZ | CCLE, pérdida de número de copias | CLL-cnd-GEN | CLL | GEN | 515,275 | 21,874 | Dirigida (bipartita) | Integrado |
| achilles_HMZ | Essentiality "bad fitness" (CRISPR/RNAi) | CLL-bfn-GEN | CLL | GEN | 49,264 | 4,822 | Dirigida (bipartita) | Integrado |
| achilles_HMZ | Essentiality "good fitness" | CLL-gfn-GEN | CLL | GEN | 49,374 | 4,822 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Increases expression | CPD-upe-GEN | CPD | GEN | 326,318 | 23,499 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Decreases expression | CPD-dwe-GEN | CPD | GEN | 311,488 | 22,173 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Affects expression (sin dirección reportada) | CPD-afe-GEN | CPD | GEN | 63,392 | 16,521 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Increases biotransformación/metabolismo | CPD-upm-GEN | CPD | GEN | 89,557 | 21,678 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Decreases biotransformación/metabolismo | CPD-dwm-GEN | CPD | GEN | 90,381 | 20,346 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Affects biotransformación/metabolismo | CPD-afm-GEN | CPD | GEN | 28,502 | 13,420 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Increases abundancia | CPD-upb-GEN | CPD | GEN | 61,463 | 13,685 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Decreases abundancia | CPD-dwb-GEN | CPD | GEN | 812 | 634 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Increases actividad | CPD-upa-GEN | CPD | GEN | 17,947 | 6,912 | Dirigida (bipartita) | Integrado |
| ctdchemgen | Decreases actividad | CPD-dwa-GEN | CPD | GEN | 14,344 | 7,077 | Dirigida (bipartita) | Integrado |
| ctdchemis | DirectEvidence = marker/mechanism ("causes") | CPD-cau-DIS | CPD | DIS | 69,051 | 9,222 | Dirigida (bipartita) | Integrado |
| ctdchemis | DirectEvidence = therapeutic ("treats") | CPD-trt-DIS | CPD | DIS | 39,515 | 9,573 | Dirigida (bipartita) | Integrado |

**Verificación de calidad ejecutada sobre las 38 filas (no estimada):**

- Filas duplicadas exactas (mismo node1, node2, tipos, relación y fuente) en las 110,690,581 aristas totales: **0**.
- Duplicados de un mismo par (node1, node2) dentro de la misma fuente+relación: **0** en las 38 combinaciones.
- "Dirección" con % bidireccional se calculó empíricamente comparando pares (A,B) vs (B,A) dentro de cada fuente+relación (se omitió en `coexpressdb` por tamaño — su simetría está garantizada por construcción en `script.py`, que usa `tuple(sorted(...))`).

---

## 2. Solapamiento semántico entre fuentes (relaciones canónicas compartidas por más de una fuente)

| Relación | Fuente A | Fuente B | Pares en A | Pares en B | Pares compartidos | % sobre A | % sobre B |
|---|---|---|---|---|---|---|---|
| GEN-ppi-GEN | omnipath | string | 84,608 | 229,073 | 23,798 | 28.1% | 10.4% |
| GEN-reg-GEN | dorothea_AB | dorothea_CD | 16,287 | 285,934 | 259 | 1.6% | 0.1% |
| GEN-upr-GEN | dorothea_AB | dorothea_CD | 6,528 | 12,610 | 3 | 0.05% | 0.02% |
| GEN-dwr-GEN | dorothea_AB | dorothea_CD | 1,097 | 3,503 | 6 | 0.5% | 0.2% |

*(Los conteos "Pares en A/B" de esta tabla difieren levemente de la Sección 1 porque aquí se cuentan pares no dirigidos únicos, es decir (A,B) y (B,A) colapsados a un mismo par, para poder comparar contra la otra fuente sin importar el sentido.)*

---

## 3. Nodos y aristas por tipo (graph_v4 completo)

| Tipo de nodo | Significado | Nodos únicos |
|---|---|---|
| GEN | Genes | 26,330 |
| DIS | Enfermedades | 23,697 |
| CPD | Compuestos/químicos | 19,629 |
| PWY | Rutas metabólicas | 2,820 |
| CLL | Líneas celulares | 1,953 |
| TIS | Tejidos | 699 |
| **Total** | | **75,128** |

| Métrica global | Valor |
|---|---|
| Bases de datos | 17 |
| Combinaciones fuente+relación | 38 |
| Tipos de relación canónicos distintos | 34 |
| Total de aristas | 110,690,581 |
| Componentes conexas | 15 (1 componente principal con 75,100 nodos = 99.96%; 14 pares aislados = 28 nodos; 0 nodos completamente solitarios) |

## 4. Esquema de tipos de nodo y relación

Ver `graph_v4_esquema.png` adjunto. agregación de las 38 combinaciones por par de tipos de nodo conectado.

---

## 5. Informe de decisiones de integración y exclusión

### Decisión 1 — GEN-ppi-GEN: omnipath vs. string (confirmada 2026-07-26)

**Hallazgo**: 84,608 pares en omnipath, 229,073 en string, solo 23,798 compartidos (28.1% / 10.4%).

**Propuesta**: mantener como fuentes separadas (no fusionar en una única relación canónica), preservando la columna `database` como trazabilidad de evidencia. El bajo solapamiento indica que ambas bases aportan evidencia mayormente complementaria. string cubre ~2.7x más interacciones que omnipath, y de las de omnipath, solo el 28% están corroboradas también en string. Fusionar perdería la distinción de qué arista viene de qué fuente sin ganar nada en deduplicación real (los duplicados exactos ya dieron 0).

**Alternativa considerada y descartada**: crear una relación `GEN-ppi-GEN` única con score combinado (como por ejemplo el máximo o promedio de confianza entre fuentes). Se descarta porque introduciría una decisión de ponderación no trivial y no está alineada con la metodología de Bioteque (que trata cada fuente como evidencia independiente).

### Decisión 2 — dorothea_AB vs. dorothea_CD (confirmada 2026-07-26)

**Hallazgo**: solapamiento casi nulo entre niveles de confianza (0.02%–1.6%).

**Propuesta**: mantener separados. AB (ChIP-seq/literatura, alta confianza) y CD (predicciones computacionales, confianza media-baja) capturan conjuntos de interacciones TF-gen mayoritariamente distintos, no la misma evidencia replicada en dos niveles. Esto ya estaba implícito en el diseño original del proyecto (documentado en runs anteriores), pero ahora queda confirmado empíricamente con datos de graph_v4, no solo por diseño.

### Decisión 3 — Corrección de documentación: `opentargets/README.txt` (aplicada 2026-07-26)

**Hallazgo**: el README interno de `opentargets` todavía decía `association_overall_direct`, pero tanto `get_data.sh` como el grafo construido usan correctamente `association_overall_indirect` (cambio aplicado deliberadamente alineado a la metodología de Bioteque). El README no se había actualizado cuando se hizo el cambio. De paso se corrigió también una segunda nota desactualizada: el README decía "Se usa https:// en wget (no ftp://)", cuando en realidad `get_data.sh` usa `ftp://` (la conexión `https://` se cuelga tras el TLS handshake en el cluster).

**Acción**: README corregido para reflejar la realidad del pipeline. No afecta a graph_v4 en sí (los datos ya estaban correctos), era puramente deuda de documentación.

### Decisión 4 — Archivos huérfanos sin uso (informativo, no bloqueante)

Se confirmó que estos archivos descargados no son referenciados por ningún `script.py` y no afectan graph_v4: `reactome/ReactomePathways.txt`, `reactome/Ensembl2Reactome_All_Levels.txt`, `cclerna_HMZ/OmicsProfiles.csv`. Quedan documentados como candidatos de limpieza de disco, sin acción obligatoria.

### Confirmaciones

- No hay aristas duplicadas exactas en graph_v4 (0 de 110,690,581).
- No hay duplicados de par dentro de ninguna fuente+relación individual.
- Las fechas de descarga de `string`, `hpa_rna_cons` y `omnipath` (enz_sub) que parecían antiguas fueron verificadas contra el servidor de origen: coinciden byte a byte con el contenido actualmente servido. No son datos obsoletos, simplemente esas fuentes no han publicado una versión más nueva.
