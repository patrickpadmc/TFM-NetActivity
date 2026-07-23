#!/usr/bin/env bash
# =============================================================================
# opentargets — get_data.sh
# Fuente: Open Targets Platform 26.03
# Contiene: GEN-ass-DIS (gene associates disease)
#
# NOTA sobre paths (post-25.03):
#   - Solo formato Parquet disponible
#   - Usar https:// en lugar de ftp:// (recomendado por la documentación oficial)
#   - cut-dirs=8 para https, ajustado desde el valor 6 del ftp original
#
# Datasets descargados:
#   - association_overall_direct: scores de asociación gen-enfermedad
#   - disease: metadata de enfermedades (nombres, ontologías EFO/MONDO)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/opentargets}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

BASE_FTP="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output"

echo "[opentargets] Descargando association_overall_direct (Parquet)..."
wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
     --continue --show-progress \
     "${BASE_FTP}/association_overall_direct"

echo "[opentargets] Descargando disease (Parquet)..."
wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
     --continue --show-progress \
     "${BASE_FTP}/disease"

echo "[opentargets] Listo. Archivos en: $OUT_DIR"
du -sh "$OUT_DIR"
