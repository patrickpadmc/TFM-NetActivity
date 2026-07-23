#!/usr/bin/env bash
# =============================================================================
# hpa_rna_cons — get_data.sh
# Fuente: Human Protein Atlas
# Archivo: rna_tissue_consensus.tsv.zip
# Contiene: GEN-upr-TIS (tissue upregulates) y GEN-dwr-TIS (downregulates)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/hpa_rna_cons}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URL="https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"

echo "[hpa_rna_cons] Descargando $URL ..."
wget --continue --show-progress -O rna_tissue_consensus.tsv.zip "$URL"

echo "[hpa_rna_cons] Descomprimiendo..."
unzip -o rna_tissue_consensus.tsv.zip
rm -f rna_tissue_consensus.tsv.zip

echo "[hpa_rna_cons] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
