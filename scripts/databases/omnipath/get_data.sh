#!/usr/bin/env bash
# =============================================================================
# omnipath — get_data.sh
# Fuente: OmniPath (archive.omnipathdb.org)
# Archivos:
#   - omnipath_webservice_interactions__latest.tsv.gz  -> ppi, dorothea
#   - omnipath_webservice_enz_sub__latest.tsv.gz       -> dephosphorylation
# Contiene: GEN-ppi-GEN, GEN-pho-GEN, GEN-dph-GEN
#
# NOTA: dorothea_AB y dorothea_CD usan el mismo archivo de interacciones.
#   Este script descarga los archivos compartidos una sola vez.
#   Los scripts de dorothea_AB y dorothea_CD hacen symlink a estos archivos.
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/omnipath}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[omnipath] Descargando interactions (ppi, phosphorylation, dorothea)..."
wget --continue --show-progress \
    -O omnipath_webservice_interactions__latest.tsv.gz \
    "https://archive.omnipathdb.org/omnipath_webservice_interactions__latest.tsv.gz"

echo "[omnipath] Descargando enzyme-substrate (dephosphorylation)..."
wget --continue --show-progress \
    -O omnipath_webservice_enz_sub__latest.tsv.gz \
    "https://archive.omnipathdb.org/omnipath_webservice_enz_sub__latest.tsv.gz"

echo "[omnipath] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
