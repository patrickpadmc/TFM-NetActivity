#!/usr/bin/env bash
# =============================================================================
# string — get_data.sh
# Fuente: STRING v12.0 (human, taxid 9606)
# NOTA (v4): human.uniprot_2_string.2018.tsv.gz fue descontinuado por STRING.
# No se descarga porque string/script.py no lo usa -- usa data/metadata/ensp2ensg.tsv.gz
# =============================================================================
set -euo pipefail
OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/string}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
echo "[string] Descargando protein.links.full..."
wget --continue --show-progress \
    -O 9606.protein.links.full.v12.0.txt.gz \
    "https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"
echo "[string] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
