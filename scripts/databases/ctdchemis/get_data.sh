#!/usr/bin/env bash
# =============================================================================
# ctdchemis — get_data.sh
# Fuente: Comparative Toxigenomics Database (CTD)
# Archivos: CTD_chemicals_diseases.tsv.gz
# Contiene: CPD-cau-DIS (compound causes disease) y CPD-trt-DIS (treats disease)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/ctdchemis}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[ctdchemis] Descargando CTD_chemicals_diseases.tsv.gz ..."
wget --continue --show-progress \
    -O CTD_chemicals_diseases.tsv.gz \
    "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz"

echo "[ctdchemis] Listo (mantenemos gzip para procesado eficiente en memoria)."
ls -lh "$OUT_DIR"
