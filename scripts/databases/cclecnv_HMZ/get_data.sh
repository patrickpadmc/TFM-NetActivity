#!/usr/bin/env bash
set -euo pipefail
OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/cclecnv_HMZ}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
DOWNLOADER="${OUT_DIR}/harmonizomedownloader.py"
wget --continue --show-progress -O "$DOWNLOADER" "https://maayanlab.cloud/Harmonizome/static/harmonizomedownloader.py"
python3 - <<'PYEOF'
import sys
sys.path.insert(0, ".")
from harmonizomedownloader import download_datasets
download_datasets(
    [('CCLE Cell Line Gene CNV Profiles', 'cclecnv')],
    ['gene_attribute_edges.txt.gz'],
)
PYEOF
mv -f "CCLE Cell Line Gene CNV Profiles/gene_attribute_edges.txt.gz" .
rm -rf "CCLE Cell Line Gene CNV Profiles"
echo "[cclecnv_HMZ] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
