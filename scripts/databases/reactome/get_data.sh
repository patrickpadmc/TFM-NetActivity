#!/usr/bin/env bash
# =============================================================================
# reactome — get_data.sh
# Fuente: Reactome v97
# Archivo: UniProt2Reactome_All_Levels.txt
# Contiene: GEN-ass-PWY (gene associates pathway)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/reactome}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URL="https://download.reactome.org/97/UniProt2Reactome_All_Levels.txt"

echo "[reactome] Descargando $URL ..."
wget --continue --show-progress -O UniProt2Reactome_All_Levels.txt "$URL"

echo "[reactome] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
