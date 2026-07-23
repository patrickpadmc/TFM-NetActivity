#!/usr/bin/env bash
# =============================================================================
# dorothea_AB — get_data.sh
# Fuente: DoRothEA via OmniPath (confidence levels A y B)
# Contiene: GEN-upr-GEN, GEN-reg-GEN, GEN-dwr-GEN
#
# Los archivos fuente son los mismos que omnipath. Este script crea
# symlinks para evitar duplicar ~GBs de datos en disco.
# Si omnipath no está descargado aún, lo descarga primero.
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/dorothea_AB}"
OMNI_DIR="$(dirname "$OUT_DIR")/omnipath"
mkdir -p "$OUT_DIR"

INTERACTIONS_FILE="${OMNI_DIR}/omnipath_webservice_interactions__latest.tsv.gz"

# Descarga omnipath si no existe
if [[ ! -f "$INTERACTIONS_FILE" ]]; then
    echo "[dorothea_AB] Archivos de omnipath no encontrados. Descargando..."
    mkdir -p "$OMNI_DIR"
    wget --continue --show-progress \
        -O "$INTERACTIONS_FILE" \
        "https://archive.omnipathdb.org/omnipath_webservice_interactions__latest.tsv.gz"
fi

# Crear symlink
ln -sf "$INTERACTIONS_FILE" "${OUT_DIR}/omnipath_webservice_interactions__latest.tsv.gz"

echo "[dorothea_AB] Symlink creado → $INTERACTIONS_FILE"
ls -lh "$OUT_DIR"
