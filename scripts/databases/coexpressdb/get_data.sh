#!/usr/bin/env bash
# =============================================================================
# coexpressdb — get_data.sh
# Fuente: COXPRESdb v22-05 (Hsa-r.c6-0)
# Contiene: GEN-cex-GEN (gene coexpresses gene)
#
# ⚠️  VERIFICAR ANTES DE EJECUTAR:
#   El archivo principal se infirió del índice del directorio. Confirma que
#   este URL es correcto antes de lanzar el job:
#   https://coxpresdb.jp/download/Hsa-r.c6-0/coex_md5/
#
# Bioteque usó (versión antigua c4-0):
#   wget --content-disposition http://coxpresdb.jp/download/Hsa-r.c4-0/coex
#   wget https://coxpresdb.jp/static/download/supportability.2014-08-19.txt
#
# Versión actual (c6-0 / v22-05):
#   Archivo de coexpresión: .zip inferido desde el .md5 que encontraste
#   Archivo de supportability: actualizado a ver8-0
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/coexpressdb}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Archivo principal de coexpresión (⚠️ confirmar URL)
COEX_URL="https://coxpresdb.jp/download/Hsa-r.c6-0/coex_md5/Hsa-r.v22-05.G16651-S235187.combat_pca.subagging.z.d.zip"
COEX_FILE="Hsa-r.v22-05.G16651-S235187.combat_pca.subagging.z.d.zip"

# Archivo de supportability (versión actualizada)
SUPP_URL="https://coxpresdb.jp/static/download/supportability.ver8-0.tar.bz2"
SUPP_FILE="supportability.ver8-0.tar.bz2"

echo "[coexpressdb] Descargando archivo de coexpresión..."
wget --continue --show-progress -O "$COEX_FILE" "$COEX_URL"

echo "[coexpressdb] Descomprimiendo coexpresión..."
unzip -o "$COEX_FILE"
rm -f "$COEX_FILE"

echo "[coexpressdb] Descargando archivo de supportability..."
wget --continue --show-progress -O "$SUPP_FILE" "$SUPP_URL"

echo "[coexpressdb] Descomprimiendo supportability..."
tar -xjf "$SUPP_FILE"
rm -f "$SUPP_FILE"

echo "[coexpressdb] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
