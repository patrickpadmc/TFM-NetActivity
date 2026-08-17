# Paquete de reproducibilidad: benchmark funcional TFM-NetActivity v1

Este directorio contiene la instantánea de resultados, particiones y controles de calidad que sustentan la evaluación final de embeddings del proyecto TFM-NetActivity. Los archivos se publican junto con `SHA256SUMS.txt`, que permite verificar su integridad.

## Estructura

### `tables/`

Resultados finales de los análisis downstream.

- `go_all_term_results.tsv`: resultados por combinación de embedding y término GO Biological Process; 616 filas de evaluación.
- `biogrid_all_fold_results.tsv`: resultados por embedding, tarea y fold externo de BioGRID; 66 filas.
- `functional_benchmark_master_summary.tsv`: resumen integrado de GO-BP y BioGRID para los 11 embeddings.
- `seccion12_tabla_completa.tsv`: evaluación estructural mediante predicción de enlaces; incluye tamaños de partición, selección de modelo, hiperparámetros, métricas de validación y métricas de test por fuente y relación.
- `seccion12_matriz_auprc.tsv`: matriz de AUPRC de test empleada en la comparación estructural.
- `seccion12_hiperparametros_umbrales.tsv`: hiperparámetros seleccionados y umbrales de decisión para la evaluación estructural.

### `splits/`

Particiones en identificadores Ensembl Gene (`ENSG`) utilizadas para prevenir fuga de información.

- `go_reference_splits_common_ensg.tsv`: particiones train/validación/holdout de 56 términos GO-BP; 82.721 ejemplos.
- `biogrid_reference_nested_splits_common_ensg.tsv`: particiones anidadas de BioGRID para letalidad sintética e interacción genética negativa; 506.085 pares.
- Los archivos `*_summary.tsv` documentan el tamaño y la composición de las particiones.
- Los archivos `*_collisions_removed.tsv` registran las colisiones ENSG eliminadas para asegurar ausencia de solapamiento entre particiones.

### `coverage/`

- `coverage_summary.json`: cobertura de los embeddings externos UniProt y GenePT sobre los genes del subgrafo común.
- `coverage_by_relation.tsv`: cobertura de embeddings por relación evaluada.

## Integridad y trazabilidad

`INVENTORY.txt` lista todos los archivos incluidos y sus tamaños.  
`SHA256SUMS.txt` contiene los checksums SHA-256 de cada archivo, para comprobar que una copia descargada coincide con la instantánea publicada.

Los scripts que generan las tablas y los benchmarks se encuentran en `scripts/evaluation/`. La definición y auditoría de `graph_v4` se documentan en `docs/graph_v4/`.

Las bases de datos fuente y los embeddings de gran tamaño no se redistribuyen en este paquete. Las fuentes originales, criterios de procesado, versiones y mapeos se describen en la documentación y en los scripts versionados del repositorio.
