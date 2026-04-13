#!/bin/bash
# Script maestro para análisis de DisGeNET y Open Targets
# Uso: bash run_analysis.sh /ruta/datos

set -e

if [ $# -lt 1 ]; then
    echo "Uso: bash run_analisis.sh <data_dir> [python_version]"
    echo ""
    echo "Ejemplo:"
    echo "  bash run_analisis.sh ~/TFM-NetActivity/data"
    echo "  bash run_analisis.sh ~/TFM-NetActivity/data python3"
    exit 1
fi

DATA_DIR="$1"
PYTHON_CMD="${2:-python3}"

# Detectar versión de Python
echo "=========================================="
echo "ANÁLISIS: DisGeNET y Open Targets"
echo "=========================================="
echo ""
echo "Verificando Python..."
$PYTHON_CMD --version || {
    echo " Python no encontrado. Intenta con:"
    echo "  bash run_analisis.sh $DATA_DIR python"
    exit 1
}

# Paths
DISGENET_FILE="$DATA_DIR/disgenet_25.4/all_gene_disease_associations.tsv.gz"
OPENTARGETS_DIR="$DATA_DIR/opentargets_26.03/association_by_datasource_direct"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""
echo "=========================================="
echo "VERIFICANDO ARCHIVOS"
echo "=========================================="

# Verificar DisGeNET
if [ -f "$DISGENET_FILE" ]; then
    echo " DisGeNET encontrado: $DISGENET_FILE"
    ANALYZE_DISGENET=1
else
    echo " DisGeNET NO encontrado en: $DISGENET_FILE"
    ANALYZE_DISGENET=0
fi

# Verificar Open Targets
if [ -d "$OPENTARGETS_DIR" ]; then
    PARQUET_COUNT=$(find "$OPENTARGETS_DIR" -name "*.parquet" 2>/dev/null | wc -l)
    echo " Open Targets encontrado: $PARQUET_COUNT archivos parquet"
    ANALYZE_OPENTARGETS=1
else
    echo " Open Targets NO encontrado en: $OPENTARGETS_DIR"
    ANALYZE_OPENTARGETS=0
fi

echo ""
echo "=========================================="
echo "EJECUTANDO ANÁLISIS"
echo "=========================================="

# Crear directorio de resultados
RESULTS_DIR="./analisis_results"
mkdir -p "$RESULTS_DIR"

# Análisis DisGeNET
if [ $ANALYZE_DISGENET -eq 1 ]; then
    echo ""
    echo "[1/2] Analizando DisGeNET..."
    $PYTHON_CMD "$SCRIPT_DIR/analisis_disgenet.py" "$DISGENET_FILE" "$RESULTS_DIR/disgenet"
    DISGENET_DONE=1
else
    echo ""
    echo "[1/2] DisGeNET: SALTADO (archivo no encontrado)"
    DISGENET_DONE=0
fi

# Análisis Open Targets
if [ $ANALYZE_OPENTARGETS -eq 1 ]; then
    echo ""
    echo "[2/2] Analizando Open Targets..."
    $PYTHON_CMD "$SCRIPT_DIR/analisis_opentargets.py" "$OPENTARGETS_DIR" "$RESULTS_DIR/opentargets"
    OPENTARGETS_DONE=1
else
    echo ""
    echo "[2/2] Open Targets: SALTADO (directorio no encontrado)"
    OPENTARGETS_DONE=0
fi

# ============================================
# RESUMEN FINAL
# ============================================
echo ""
echo "=========================================="
echo "ANÁLISIS COMPLETADO"
echo "=========================================="

echo ""
echo "📊 Resultados guardados en: $RESULTS_DIR/"
echo ""

if [ $DISGENET_DONE -eq 1 ]; then
    echo "DisGeNET (gene-disease-score):"
    echo "  - descriptives.csv"
    echo "  - gene_disease_score.csv"
    echo "  - gene_disease_score_top1000.csv"
    echo "  - score_distribution.csv"
    echo "  - sources_distribution.csv"
    echo ""
fi

if [ $OPENTARGETS_DONE -eq 1 ]; then
    echo "Open Targets (target-disease-score):"
    echo "  - descriptives.csv"
    echo "  - target_disease_score.csv"
    echo "  - target_disease_score_top1000.csv"
    echo "  - score_distribution.csv"
    echo ""
fi

# Crear archivo de resumen
cat > "$RESULTS_DIR/SUMMARY.txt" << EOF
ANÁLISIS DE BASES DE DATOS BIOINFORMÁTICAS
===========================================

Fecha: $(date)

DisGeNET 25.4
- Tipo: Gene-Disease Associations
- Archivo: $DISGENET_FILE
- Estado: $([ $DISGENET_DONE -eq 1 ] && echo "COMPLETADO" || echo "NO ANALIZADO")

Open Targets Platform 26.03
- Tipo: Target-Disease Evidence (Parquet format)
- Directorio: $OPENTARGETS_DIR
- Estado: $([ $OPENTARGETS_DONE -eq 1 ] && echo "COMPLETADO" || echo "NO ANALIZADO")

Próximos pasos:
1. Revisar los archivos CSV generados
2. Comparar gene-disease scores vs target-disease scores
3. Usar los archivos top1000 para análisis inicial
4. Explorar distribuciones de scores en ambas bases

EOF

echo "✓ Resumen guardado en: $RESULTS_DIR/SUMMARY.txt"
echo ""
echo "=========================================="
echo "¡LISTO PARA ANÁLISIS POSTERIOR!"
echo "=========================================="
echo ""