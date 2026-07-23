#!/usr/bin/env bash
# =============================================================================
# dorothea_CD — get_data.sh
# Fuente: DoRothEA via OmniPath (confidence levels C y D)
# Contiene: GEN-upr-GEN, GEN-reg-GEN, GEN-dwr-GEN
#
# Igual que dorothea_AB pero para niveles de confianza C y D.
# El filtrado por nivel se hace en script.py.
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/dorothea_CD}"
OMNI_DIR="$(dirname "$OUT_DIR")/omnipath"
mkdir -p "$OUT_DIR"

INTERACTIONS_FILE="${OMNI_DIR}/omnipath_webservice_interactions__latest.tsv.gz"

if [[ ! -f "$INTERACTIONS_FILE" ]]; then
    echo "[dorothea_CD] Archivos de omnipath no encontrados. Descargando..."
    mkdir -p "$OMNI_DIR"
    wget --continue --show-progress \
        -O "$INTERACTIONS_FILE" \
        "https://archive.omnipathdb.org/omnipath_webservice_interactions__latest.tsv.gz"
fi

ln -sf "$INTERACTIONS_FILE" "${OUT_DIR}/omnipath_webservice_interactions__latest.tsv.gz"

echo "[dorothea_CD] Symlink creado → $INTERACTIONS_FILE"
ls -lh "$OUT_DIR"
