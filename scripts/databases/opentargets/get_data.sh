#!/usr/bin/env bash
# =============================================================================
# opentargets — get_data.sh
# Fuente: Open Targets Platform 26.06
# Contiene: GEN-ass-DIS (gene associates disease) -- version INDIRECT (v4)
#
# NOTA (v4): HTTPS a ftp.ebi.ac.uk se cuelga tras el TLS handshake en este
# cluster (nunca responde HTTP). FTP puro funciona bien -- usamos ftp://.
#
# Datasets descargados:
#   - association_overall_indirect: scores gen-enfermedad (incl. inferidas via EFO/MONDO)
#   - disease: metadata de enfermedades (nombres, ontologias EFO/MONDO)
# =============================================================================
set -euo pipefail
OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/opentargets}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
BASE_FTP="ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output"

echo "[opentargets] Descargando association_overall_indirect (Parquet, ftp)..."
wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
     --continue --show-progress \
     "${BASE_FTP}/association_overall_indirect/"

echo "[opentargets] Descargando disease (Parquet, ftp)..."
wget --recursive --no-parent --no-host-directories --cut-dirs 8 \
     --continue --show-progress \
     "${BASE_FTP}/disease/"

echo "[opentargets] Listo. Archivos en: $OUT_DIR"
du -sh "$OUT_DIR"
find "$OUT_DIR" -name "*.parquet" | wc -l
