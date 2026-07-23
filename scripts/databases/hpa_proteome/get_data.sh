#!/usr/bin/env bash
# =============================================================================
# hpa_proteome — get_data.sh
# Fuente: Human Protein Atlas
# Archivo: normal_ihc_data.tsv.zip
# Contiene: GEN-pdf-TIS (protein deficiency) y GEN-pab-TIS (protein abundance)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/hpa_proteome}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URL="https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"

echo "[hpa_proteome] Descargando $URL ..."
wget --continue --show-progress -O normal_ihc_data.tsv.zip "$URL"

echo "[hpa_proteome] Descomprimiendo..."
unzip -o normal_ihc_data.tsv.zip
rm -f normal_ihc_data.tsv.zip

echo "[hpa_proteome] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
