#!/usr/bin/env bash
# =============================================================================
# ctddisease — get_data.sh
# Fuente: Comparative Toxigenomics Database (CTD)
# Archivos: CTD_genes_diseases.tsv.gz
# Contiene: DIS-ass-GEN (gene associates disease)
#
# NOTA: ctdchemis también usa CTD_chemicals_diseases.tsv.gz.
# Ambas databases comparten la misma fuente; si ya descargaste ctdchemis
# puedes copiar los archivos en vez de re-descargarlos.
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/ctddisease}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[ctddisease] Descargando CTD_genes_diseases.tsv.gz ..."
wget --continue --show-progress \
    -O CTD_genes_diseases.tsv.gz \
    "https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz"

echo "[ctddisease] Listo (mantenemos gzip para procesado eficiente en memoria)."
ls -lh "$OUT_DIR"
