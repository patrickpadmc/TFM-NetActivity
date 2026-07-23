#!/usr/bin/env bash
# =============================================================================
# repohub — get_data.sh
# Fuente: The Drug Repurposing Hub (Broad Institute)
# Archivo: repo-drug-annotation-20200324.txt
# Contiene: CPD-int-GEN (compound interacts gene)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/repohub}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URL="https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20200324.txt"

echo "[repohub] Descargando $URL ..."
wget --continue --show-progress -O repo-drug-annotation-20200324.txt "$URL"

echo "[repohub] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
