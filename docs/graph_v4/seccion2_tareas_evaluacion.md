# Sección 2 — Selección de tareas de evaluación

**Fecha**: 2026-07-26
**Basado en**: graph_v4 congelado (commit `699d008`), tabla de auditoría de Sección 1.

Todos los conteos de esta sección provienen directamente de la tabla de auditoría de la Sección 1 (verificada contra `graph_edges.tsv`, 110,690,581 aristas). Ningún número es estimado.

---

## 1. Umbral de soporte

El protocolo sugiere un mínimo de 500 aristas positivas **después de la división de datos**. El proyecto ya tiene establecida una convención de split 90/10 (`seed42_test10`, documentada en ejecuciones previas del proyecto). Bajo ese split, el conjunto de test recibe ~10% de las aristas de cada relación.

Para que el test tenga ≥500 positivos: **total de aristas de la relación ≥ 500 / 0.10 = 5,000**.

Se adopta **5,000 aristas totales** como umbral mínimo de soporte para que una relación sea evaluable formalmente. Dado que la escala real del grafo va de 749 hasta 105,907,896 aristas por relación, este umbral es deliberadamente bajo. Solo descarta los casos extremos de escasez real; no pretende ser el criterio para elegir tareas *primarias*.

Cuando una relación tiene evidencia de más de una fuente, el conteo de aristas usado es el de **pares únicos tras deduplicar entre fuentes** (no la suma cruda de filas, que contaría dos veces la misma interacción biológica si dos fuentes la reportan).

---

## 2. Tabla de relaciones: aristas y decisión de uso

*(Familia asignada según el par de tipos de nodo. "Aristas (dedupe)" = aristas únicas tras eliminar duplicados entre fuentes, cuando aplica.)*

| Familia | Relación | Fuente(s) | Aristas (dedupe) | ≥5,000 | Decisión |
|---|---|---|---|---|---|
| Gen-enfermedad | GEN-ass-DIS | opentargets | 283,977 | Sí | **Primaria** |
| Gen-tejido | GEN-ass-TIS | jensentissuecurated | 315,368 | Sí | **Primaria** |
| Gen-tejido | GEN-pab-TIS | hpa_proteome | 77,568 | Sí | Secundaria |
| Gen-tejido | GEN-pdf-TIS | hpa_proteome | 103,702 | Sí | Secundaria |
| Gen-tejido | GEN-upr-TIS | hpa_rna_cons | 12,750 | Sí | Secundaria |
| Gen-tejido | GEN-dwr-TIS | hpa_rna_cons | 12,750 | Sí | Secundaria |
| Gen-vía | GEN-ass-PWY | reactome | 152,359 | Sí | **Primaria** |
| Compuesto-enfermedad | CPD-trt-DIS | ctdchemis | 39,515 | Sí | **Primaria** |
| Compuesto-enfermedad | CPD-cau-DIS | ctdchemis | 69,051 | Sí | Secundaria |
| Compuesto-gen | CPD-int-GEN | repohub | 10,279 | Sí | **Primaria** |
| Compuesto-gen | CPD-upe-GEN | ctdchemgen | 326,318 | Sí | Secundaria |
| Compuesto-gen | CPD-dwe-GEN | ctdchemgen | 311,488 | Sí | Secundaria |
| Compuesto-gen | CPD-afe-GEN | ctdchemgen | 63,392 | Sí | Secundaria |
| Compuesto-gen | CPD-upm-GEN | ctdchemgen | 89,557 | Sí | Secundaria |
| Compuesto-gen | CPD-dwm-GEN | ctdchemgen | 90,381 | Sí | Secundaria |
| Compuesto-gen | CPD-afm-GEN | ctdchemgen | 28,502 | Sí | Secundaria |
| Compuesto-gen | CPD-upb-GEN | ctdchemgen | 61,463 | Sí | Secundaria |
| Compuesto-gen | CPD-dwb-GEN | ctdchemgen | 812 | **No** | **Excluida** (soporte insuficiente) |
| Compuesto-gen | CPD-upa-GEN | ctdchemgen | 17,947 | Sí | Secundaria |
| Compuesto-gen | CPD-dwa-GEN | ctdchemgen | 14,344 | Sí | Secundaria |
| Gen-gen | GEN-ppi-GEN | omnipath + string | 289,883 | Sí | **Primaria** *(ver Sección 4) |
| Gen-gen | GEN-cex-GEN | coexpressdb | 105,907,896 | Sí | Secundaria (ver nota) |
| Gen-gen | GEN-pho-GEN | omnipath | 27,639 | Sí | Secundaria |
| Gen-gen | GEN-dph-GEN | omnipath | 749 | **No** | **Excluida** (soporte insuficiente) |
| Gen-gen | GEN-reg-GEN | dorothea_AB + dorothea_CD | 301,962 | Sí | Secundaria * (ver Sección 4) |
| Gen-gen | GEN-upr-GEN | dorothea_AB + dorothea_CD | 19,135 | Sí | Secundaria * (ver Sección 4) |
| Gen-gen | GEN-dwr-GEN | dorothea_AB + dorothea_CD | 4,594 | **No** | **Excluida** (soporte insuficiente) |
| Línea celular-gen | CLL-mut-GEN | cclemut_HMZ | 107,951 | Sí | **Primaria** |
| Línea celular-gen | CLL-upr-GEN | cclerna_HMZ | 401,250 | Sí | Secundaria |
| Línea celular-gen | CLL-dwr-GEN | cclerna_HMZ | 401,250 | Sí | Secundaria |
| Línea celular-gen | CLL-cnu-GEN | cclecnv_HMZ | 508,053 | Sí | Secundaria |
| Línea celular-gen | CLL-cnd-GEN | cclecnv_HMZ | 515,275 | Sí | Secundaria |
| Línea celular-gen | CLL-bfn-GEN | achilles_HMZ | 49,264 | Sí | Secundaria |
| Línea celular-gen | CLL-gfn-GEN | achilles_HMZ | 49,374 | Sí | Secundaria |

**Nota sobre GEN-cex-GEN**: clasificada como secundaria pese a superar ampliamente el umbral (105.9M aristas, 95.7% del grafo). Motivo: no es soporte lo que falta, es (a) que la coexpresión es evidencia estadística indirecta (correlación transcripcional), no interacción física curada. Semánticamente más débil que GEN-ppi-GEN para el propósito de "predicción de enlaces biológicamente interpretable"; y (b) evaluarla a esa escala exige decisiones de muestreo de negativos y cómputo muy distintas al resto de tareas, lo que la hace poco comparable en igualdad de condiciones frente a las demás familias. Se mantiene íntegra en el grafo para el entrenamiento de embeddings (es la columna vertebral estructural del grafo), pero no se usa como tarea de comparación formal.

**Resumen de exclusiones**: 3 relaciones excluidas por soporte insuficiente (< 5,000 aristas): `GEN-dph-GEN` (749), `CPD-dwb-GEN` (812), `GEN-dwr-GEN` combinada dorothea_AB+CD (4,594). Ninguna se elimina del grafo, solo quedan fuera de la evaluación formal de predicción de enlaces.

---

## 3. Tareas primarias propuestas

Siete tareas, una por familia biológica, cubriendo todo el espectro de tipos de nodo del grafo (GEN, DIS, TIS, PWY, CPD, CLL):

| # | Familia | Relación | Aristas | Justificación biológica |
|---|---|---|---|---|
| 1 | Gen-enfermedad | **GEN-ass-DIS** | 283,977 | Única relación de la familia; es la tarea más citada en literatura de embeddings biomédicos (descubrimiento de dianas terapéuticas). Fuente única (opentargets), sin riesgo de fuga entre fuentes. |
| 2 | Gen-tejido | **GEN-ass-TIS** | 315,368 | De las 5 candidatas, es la única basada en curación directa de literatura (Jensen Lab), no en un umbral estadístico arbitrario (a diferencia de "top/bottom 250 genes por z-score" de hpa_rna_cons). Mayor volumen de la familia. Fuente única. |
| 3 | Gen-vía | **GEN-ass-PWY** | 152,359 | Única relación de la familia. Pertenencia a rutas metabólicas es un dato estructuralmente central en biología de sistemas. Fuente única. |
| 4 | Compuesto-enfermedad | **CPD-trt-DIS** | 39,515 | Confirmada 2026-07-26. Se prefiere sobre CPD-cau-DIS (69,051 aristas, mayor volumen) porque "treats" es la relación clínicamente accionable. Predecir nuevas relaciones fármaco-trata-enfermedad es la tarea central del reposicionamiento de fármacos (drug repurposing), la aplicación más citada de embeddings de grafos biomédicos. CPD-cau-DIS (toxicidad/causalidad) es una tarea válida pero de motivación distinta (farmacovigilancia). |
| 5 | Compuesto-gen | **CPD-int-GEN** | 10,279 | Confirmada 2026-07-26. Es la interacción fármaco-diana más directa y semánticamente limpia de las 11 candidatas (interacción reportada directamente, no inferida de un umbral de expresión/metilación/actividad). Corresponde al benchmark clásico de "Drug-Target Interaction (DTI) prediction" en la literatura de ML. Su volumen (10,279) es el menor de las 7 primarias propuestas, pero más que duplica el umbral de soporte. |
| 6 | Gen-gen | **GEN-ppi-GEN** | 289,883 (dedupe) | Confirmada 2026-07-26. Interacción proteína-proteína física, el benchmark de predicción de enlaces más establecido en bioinformática de redes (comparable a evaluaciones estándar sobre STRING/BioGRID). Preferida sobre GEN-cex-GEN por ser evidencia más directa/interpretable (ver nota de exclusión de cex arriba) y computacionalmente tratable a esta escala. **Requiere el procedimiento de deduplicación entre fuentes descrito en la Sección 4 antes de dividir en train/test.** |
| 7 | Línea celular-gen | **CLL-mut-GEN** | 107,951 | Confirmada 2026-07-26. Señal binaria limpia (mutación Hotspot/Damaging presente/ausente), sin discretización de un umbral continuo (a diferencia de las variantes up/down de expresión y CNV de la misma familia). Relevante para oncología de precisión (predicción de qué genes están mutados en una línea celular). |

Las 7 tareas primarias quedan confirmadas. Las tareas 1, 2 y 3 no tenían alternativa real dentro de su familia (relación única); las tareas 4, 5, 6 y 7 se eligieron entre varias candidatas viables de su misma familia.

---

## 4. Riesgos de fuga de información entre fuentes equivalentes

### 4.1 Principio general (ya establecido en el proyecto)

Cualquier embedding usado para evaluación debe calcularse sobre un grafo que excluya **todas** las aristas de test de **todas** las relaciones evaluadas, no solo la relación puntual que se está midiendo en cada momento. Esto ya estaba documentado como principio del proyecto ("Data leakage: los embeddings no pueden evaluarse sobre aristas retenidas si se calcularon sobre el grafo completo") y se reafirma aquí como base para la Sección 3.

### 4.2 Riesgo específico detectado: relaciones con más de una fuente

Tres relaciones canónicas de graph_v4 se alimentan de más de una base de datos, y la Sección 1 ya midió el solapamiento real entre ellas:

| Relación | Fuentes | Pares compartidos | % del total |
|---|---|---|---|
| GEN-ppi-GEN | omnipath, string | 23,798 | 8.2% de los 289,883 pares únicos |
| GEN-reg-GEN | dorothea_AB, dorothea_CD | 259 | 0.09% de los 301,962 pares únicos |
| GEN-upr-GEN | dorothea_AB, dorothea_CD | 3 | 0.02% de los 19,135 pares únicos |

**El riesgo concreto**: si el split de train/test se hace fila por fila sobre `graph_edges.tsv` (que mantiene omnipath y string como filas separadas para el mismo par gen-gen), es posible que la fila de `omnipath` para el par (A,B) caiga en test mientras la fila de `string` para ese mismo par (A,B) —la misma interacción biológica, solo que registrada por otra fuente— quede en train. El modelo vería literalmente esa arista durante el entrenamiento (con otro `database` en la columna, pero el mismo par de nodos) y "prediría" correctamente el enlace retenido por memorización directa, no por generalización real. Esto inflaría artificialmente el AUPRC/AUROC de esa tarea.

**Mitigación obligatoria antes de dividir GEN-ppi-GEN** (y, si se usan como secundarias, GEN-reg-GEN/GEN-upr-GEN): colapsar a pares (node1, node2) únicos **entre todas las fuentes** antes del split, no después. Es decir, deduplicar primero, dividir en train/test después. El grafo de entrenamiento de embeddings puede seguir usando las filas separadas por fuente para todo lo demás (preserva trazabilidad, como se decidió en la Sección 1); la deduplicación pre-split solo aplica al conjunto de aristas positivas de esa tarea de evaluación específica.

### 4.3 Relaciones sin este riesgo

Las 4 tareas primarias restantes (GEN-ass-DIS, GEN-ass-TIS, GEN-ass-PWY, CPD-trt-DIS) y CPD-int-GEN, CLL-mut-GEN provienen de una única fuente cada una, no existe una segunda base de datos aportando la misma relación canónica, por lo que no hay riesgo de fuga cruzada entre fuentes en estas seis tareas. Sí aplica igualmente el principio general de la Sección 4.1 (excluir aristas de test del grafo usado para generar embeddings).

### 4.4 Validación cruzada entre fuentes (mencionada por completitud, no forma parte del alcance actual)

El solapamiento medido en la Sección 1 (28.1%/10.4% para PPI) podría reutilizarse como diseño explícito de validación cruzada entre fuentes (por ejemplo: entrenar solo con string, evaluar cuántos pares de omnipath se recuperan, como proxy de "evidencia independiente"). Esto sería una pregunta de investigación distinta y válida, pero no es la evaluación de predicción de enlaces estándar que cubre esta sección, se deja mencionada como posible trabajo futuro, no como parte del protocolo actual.
