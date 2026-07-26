# Sección 3 — Particiones, negativos y protocolo común

## 1. Objetivo

Definir, para las 7 tareas primarias seleccionadas en la Sección 2, particiones train/val/test reproducibles y sin fuga de información, generar negativos válidos para cada partición, y construir un único grafo de entrenamiento combinado sobre el que se calcularán los embeddings (Spectral Embedding, Node2Vec, HIN2Vec, GAE, GAT-AE), garantizando que ningún nodo quede desconectado de ese grafo por efecto de las particiones.

## 2. Decisiones de diseño confirmadas

| # | Decisión | Confirmado |
|---|---|---|
| 1 | Grafo de entrenamiento para embeddings: **único grafo combinado** (excluye simultáneamente val+test de las 7 tareas), no un grafo separado por tarea | 2026-07-26 |
| 2 | Ratio de negativos: **1:1 en train y val** (aprendizaje balanceado), **1:10 en test** (evaluación realista, evita inflar AUPRC) | 2026-07-26 |
| 3 | Split estricto por nodos (cold-start, sin solapamiento de genes entre train/val/test): **solo GEN-ass-DIS** | 2026-07-26 |
| 4 | Split estándar (por aristas) para las 6 tareas restantes: **70% train / 10% val / 20% test** | 2026-07-26 |

## 3. Metodología

### 3.1 Deduplicación y canonicalización

Para cada tarea se extraen los pares únicos (node1, node2) de sus relaciones asociadas. En `GEN-ppi-GEN` (única relación simétrica entre las 7 tareas primarias) los pares se canonicalizan ordenando lexicográficamente antes de deduplicar, evitando contar `(A,B)` y `(B,A)` como aristas distintas, esto redujo GEN-ppi-GEN de 313,681 filas crudas a 289,883 pares únicos. El resto de tareas son bipartitas (tipos de nodo distintos en cada extremo) y no requieren canonicalización.

### 3.2 Preservación de conectividad (dos pasos)

**Paso 1: protección por grado global.** Se calculó el grado de cada uno de los 75,128 nodos de graph_v4 contando sus apariciones en **las 110,690,581 aristas del grafo completo** (no solo dentro de la tarea evaluada). 16,908 nodos (22.5% del total) tienen grado global exactamente 1: su única conexión en todo el grafo. Cualquier arista positiva con un extremo en esa situación queda forzada a train en su tarea, nunca es candidata a val/test.

**Paso 2: rescate de nodos aislados residuales.** Tras aplicar el paso 1 y excluir del grafo combinado las particiones val+test de las 7 tareas, algunos nodos con grado global ≥2 pueden quedar en grado 0 si *todas* sus aristas cayeron en held-out de una o más tareas simultáneamente. Se detectan estos nodos comparando el conjunto de nodos del grafo de entrenamiento combinado contra el conjunto completo de 75,128 nodos, y para cada nodo aislado se recupera una de sus aristas retenidas (de cualquiera de las 7 tareas) de vuelta a train.

En la ejecución real (job 308367): el paso 1 dejó **679 nodos** en grado 0 tras construir el grafo combinado; el paso 2 los rescató a los 679, dejando **0 nodos aislados** en el grafo final de entrenamiento.

### 3.3 Muestreo de negativos

Rejection sampling vectorizado por lotes: se extraen candidatos aleatorios del pool de nodos compatible por tipo (p. ej. GEN × DIS para GEN-ass-DIS), se descartan los que coinciden con un positivo real (de cualquier split) o con un negativo ya usado en otro split de la misma tarea, y se repite hasta completar la cuota. Esto garantiza que ningún negativo de train/val/test sea en realidad un enlace verdadero, y que no haya negativos duplicados entre splits.

### 3.4 Split estricto (cold-start) — solo GEN-ass-DIS

En vez de dividir aristas, se dividen los **genes**: cada gen (y todas sus aristas gen–enfermedad) se asigna íntegramente a train, val o test, de forma que ningún gen de test haya sido visto en train. 13 genes con grado global 1 quedan protegidos en train. El resto se reparte 70/10/20 por conteo de genes.

## 4. Resultados

### 4.1 Splits estándar (por arista), 7 tareas primarias

Fuente: `data/processed/evaluation/graph_v4/splits_summary.tsv` (valores finales, tras el paso de rescate).

| Tarea | Positivos totales | Protegidos (grado 1) | Train | Val | Test | Neg. train (1:1) | Neg. val (1:1) | Neg. test (1:10) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GEN-ass-DIS | 283,977 | 9,285 | 199,171 | 28,181 | 56,625 | 198,784 | 28,398 | 567,950 |
| GEN-ass-TIS | 315,368 | 85 | 220,765 | 31,530 | 63,073 | 220,757 | 31,537 | 630,740 |
| GEN-ass-PWY | 152,359 | 287 | 106,676 | 15,220 | 30,463 | 106,651 | 15,236 | 304,720 |
| CPD-trt-DIS | 39,515 | 1,650 | 27,720 | 3,920 | 7,875 | 27,660 | 3,952 | 79,030 |
| CPD-int-GEN | 10,279 | 2,288 | 7,382 | 921 | 1,976 | 7,195 | 1,028 | 20,560 |
| GEN-ppi-GEN | 289,883 | 23 | 202,919 | 28,987 | 57,977 | 202,918 | 28,988 | 579,770 |
| CLL-mut-GEN | 107,951 | 646 | 75,577 | 10,791 | 21,583 | 75,566 | 10,795 | 215,900 |

Todas las tareas alcanzan el 70/10/20 objetivo con desviaciones menores al 0.1 punto porcentual, absorbidas por el redondeo y por las aristas protegidas (que solo pueden desplazar masa hacia train, nunca hacia val/test).

### 4.2 Split estricto (cold-start), GEN-ass-DIS

| | Genes | Aristas | Negativos |
|---|---:|---:|---:|
| Train | 7,877 (protegidos: 13) | 194,783 | 194,783 (1:1) |
| Val | 1,125 | 29,252 | 29,252 (1:1) |
| Test | 2,251 | 59,942 | 599,420 (1:10) |

Archivos en `data/processed/evaluation/graph_v4/GEN-ass-DIS/{train,val,test}_strict.tsv` y `negatives_{train,val,test}_strict.tsv`.

### 4.3 Grafo de entrenamiento combinado (para cálculo de embeddings)

| Métrica | Valor |
|---|---:|
| Aristas totales en graph_v4 | 110,690,581 |
| Aristas excluidas (val+test de las 7 tareas, tras rescate) | 366,186 |
| Aristas en el grafo de entrenamiento combinado | 110,324,395 (99.67% del total) |
| Aristas rescatadas (paso 2) | 679 |
| Nodos en el grafo de entrenamiento | 75,128 / 75,128 (100%) |
| Componentes conexas | 23 |
| Componente principal | 75,084 nodos (99.94%) |
| Nodos completamente aislados (grado 0) | 0 |

El aumento de 15 a 23 componentes respecto al grafo completo (Sección 1: 15 componentes, 28 nodos en 14 pares aislados) es la consecuencia esperada de retirar 366,186 aristas para evaluación: al quitar aristas, la conectividad solo puede mantenerse o empeorar, nunca mejorar. Ningún nodo queda sin ninguna arista (0 nodos en grado 0), que era la condición mínima exigida para que todos los métodos de embedding (incluidos los que requieren que cada nodo tenga al menos una conexión, como Node2Vec) puedan generar un vector para los 75,128 nodos.

Archivo: `data/processed/evaluation/graph_v4/graph_edges_train_only.tsv`.

### 4.4 Verificación anti-fuga

Para cada una de las 7 tareas se comprobaron 4 controles, y un 5º control global:

1. Sin solapamiento de pares positivos entre train/val/test.
2. Ningún negativo coincide con un positivo real (de cualquier split).
3. Sin negativos duplicados entre splits.
4. Ninguna arista protegida (grado global 1) terminó en val/test.
5. (Global) Cero nodos con grado 0 y cero componentes de tamaño 1 en el grafo de entrenamiento combinado.

**Resultado: los 29 controles (4 × 7 tareas + 1 global) pasaron → `VERIFICACION GLOBAL: TODO OK`.**

Detalle completo en `data/processed/evaluation/graph_v4/leakage_report.txt`.

## 5. Código de evaluación común

`scripts/evaluation/common_eval.py`, compartido por los 5 métodos de embedding para que la comparación entre ellos no dependa de diferencias de implementación en la evaluación:

- `build_features()`: producto de Hadamard + diferencia absoluta entre los embeddings de los dos nodos de un par, concatenados (mínimo exigido por el protocolo).
- `train_classifier()`: regresión logística, con selección de hiperparámetro C por AUPRC en validación.
- `evaluate()`: AUPRC (métrica primaria) y AUROC (métrica secundaria) en test.
- `ranked_evaluation()`: Hits@k y MRR, rankeando cada positivo de test contra un pool de negativos muestreados (no contra el universo completo de candidatos, inviable a esta escala, hasta 26,330 × 23,697 pares posibles por tarea).

**Decisión pendiente de validar**: el tamaño del pool de negativos para `ranked_evaluation()` está en 99 por defecto (protocolo estándar "1 positivo vs N negativos" de la literatura de link prediction). Afecta la escala absoluta de Hits@k, no a AUPRC/AUROC. A confirmar antes de usarse en resultados finales del TFM, cuando se implemente la evaluación de ranking (Hits@10/MRR son métricas secundarias).

## 6. Tabla de resultados (vacía)

`results_table_v4.tsv`: 40 filas (8 combinaciones tarea×split_type × 5 métodos de embedding), columnas `task, split_type, embedding_method, feature_type, auprc, auroc, hits_at_10, mrr, notes`, métricas vacías — a rellenar en las secciones siguientes del protocolo conforme se calculen los embeddings y se evalúen.

## 7. Archivos generados

- `scripts/evaluation/build_splits_v4.py` — script de generación de particiones/negativos/grafo combinado.
- `scripts/evaluation/job_build_splits_v4.sh` — job SLURM (partición `short`, 128G, 8 cpus, 3h).
- `scripts/evaluation/common_eval.py` — código de evaluación común.
- `data/processed/evaluation/graph_v4/{TASK}/{train,val,test}.tsv` — particiones por tarea (7 tareas).
- `data/processed/evaluation/graph_v4/{TASK}/negatives_{train,val,test}.tsv` — negativos por tarea.
- `data/processed/evaluation/graph_v4/GEN-ass-DIS/{train,val,test}_strict.tsv` y `negatives_*_strict.tsv` — split estricto cold-start.
- `data/processed/evaluation/graph_v4/graph_edges_train_only.tsv` — grafo combinado de entrenamiento (110,324,395 aristas, 75,128 nodos).
- `data/processed/evaluation/graph_v4/splits_summary.tsv`, `manifest_splits.json`, `leakage_report.txt`.
- `results_table_v4.tsv` — tabla de resultados vacía.

## 8. Notas y limitaciones

- El split estricto por nodos de GEN-ass-DIS reduce el número de aristas de test evaluables (59,942) porque agrupa por gen en vez de por arista; es intencional (mide generalización a genes nunca vistos, no memorización de aristas).
- Los 28 nodos aislados del componente principal identificados en la Sección 1 (14 pares, escasez genuina de datos) permanecen en el grafo de entrenamiento combinado con su única arista intacta (no participan en ninguna de las 7 tareas evaluadas, así que no fueron afectados por las particiones).
