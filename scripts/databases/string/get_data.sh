#!/usr/bin/env bash
# =============================================================================
# string — get_data.sh
# Fuente: STRING v12.0 (human, taxid 9606)
# Archivos:
#   - 9606.protein.links.full.v12.0.txt.gz   -> interacciones PPI con scores
#   - human.uniprot_2_string.2018.tsv.gz      -> mapeo UniProt -> STRING IDs
# Contiene: GEN-ppi-GEN
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/string}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[string] Descargando protein.links.full..."
wget --continue --show-progress \
    -O 9606.protein.links.full.v12.0.txt.gz \
    "https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"

echo "[string] Descargando mapeo UniProt → STRING..."
wget --continue --show-progress \
    -O human.uniprot_2_string.2018.tsv.gz \
    "https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz"

echo "[string] Listo (archivos mantenidos en gzip). Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
