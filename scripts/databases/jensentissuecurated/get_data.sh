#!/usr/bin/env bash
# =============================================================================
# jensentissuecurated — get_data.sh
# Fuente: Jensen Lab (COMPARTMENTS)
# Archivo: human_tissue_knowledge_full.tsv
# Contiene: GEN-ass-TIS (gene associates tissue)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/jensentissuecurated}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URL="https://download.jensenlab.org/human_tissue_knowledge_full.tsv"

echo "[jensentissuecurated] Descargando $URL ..."
wget --continue --show-progress -O human_tissue_knowledge_full.tsv "$URL"

echo "[jensentissuecurated] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
