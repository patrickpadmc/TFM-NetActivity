#!/bin/bash
# # Script para descargar DisGeNET 25.4 y Open Targets 26.03
# Uso: bash download_databases.sh /ruta/donde/guardar/datos
# Ej: bash download_databases.sh ~/TFM-NetActivity/data

set -e # Salir si hay error

if [ $# -eq 0 ]; then
    echo "Error: debes especificar la ruta de destino"
    echo "Uso: bash download_databases.sh /ruta/destino"
    exit 1
fi

DEST_DIR="$1"
mkdir -p "$DEST_DIR"

echo "======================================"
echo "Descargando bases de datos al cluster"
echo "Directorio destino: $DEST_DIR"
echo "======================================"

# ======================================
# DisGeNET 25.4
# ======================================

echo "[1/2] Descargando DisGeNET 25.4 (gene-disease associations)..."
DISGENET_DIR="$DEST_DIR/disgenet_25.4"
mkdir -p "$DISGENET_DIR"

DISGENET_URL="https://www.disgenet.com/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"

cd "$DISGENET_DIR"
echo "Descargando all_gene_disease_associations.tsv.gz..."
wget -q --show-progress "$DISGENET_URL" || {
    echo " Descarga directa falló. Asegúrate de:"
    echo " 1. Tener conexión a internet."
    echo " 2. Estar registrado en disgenet.com (para versión full)."
    echo " 3. Si necesitas versión completa, descargarla desde:"
    echo " https://www.disgenet.com/downloads"
    exit 1
}
echo " DisGeNET 25.4 descargado "

# ======================================
# Open Targets Platform 26.03
# ======================================
echo "[2/2] Descargando Open Targets 26.03 (gene-disease associations)..."
OT_DIR="$DEST_DIR/opentargets_26.03"
mkdir -p "$OT_DIR"

# Descargar las asociaciones principales (associationByOverallIndirect)
# Ahora están en formato parquet por datasource (antes JSON)
# Usaremos la ruta principal de asociaciones

cd "$OT_DIR"

echo " Descargando associations split by datasource..."
FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output"

# Crear subdirectorio
mkdir -p "association_by_datasource_direct"
cd "association_by_datasource_direct"

# Usar lftp si está disponible, si no usar wget
if command -v lftp &> /dev/null; then
    echo "  → Usando lftp para descarga eficiente..."
    lftp -e "mirror --continue --parallel=4 $FTP_BASE/association_by_datasource_direct . && quit" 2>/dev/null || {
        echo "    lftp fallo, cambiando a wget..."
        wget -r -q --show-progress --no-parent --no-host-directories \
            --cut-dirs=8 "$FTP_BASE/association_by_datasource_direct/"
    }
else
    echo "  → Descargando con wget (más lento)..."
    # wget recursivo - descarga todos los parquets
    wget -r -q --show-progress --no-parent --no-host-directories \
        --cut-dirs=8 "$FTP_BASE/association_by_datasource_direct/" || {
        echo "  ⚠️  Descarga parcial. Verifica la conexión."
    }
fi

echo "  Open Targets 26.03 descargado "
 
echo ""
echo "=========================================="
echo "Descarga completada"
echo "=========================================="
echo ""
echo "Ubicaciones de los datos:"
echo "  DisGeNET:    $DISGENET_DIR"
echo "  Open Targets: $OT_DIR"
echo ""
