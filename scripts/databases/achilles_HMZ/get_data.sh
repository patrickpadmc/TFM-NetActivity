#!/usr/bin/env bash
set -euo pipefail
OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/achilles_HMZ}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
DOWNLOADER="${OUT_DIR}/harmonizomedownloader.py"
wget --continue --show-progress -O "$DOWNLOADER" "https://maayanlab.cloud/Harmonizome/static/harmonizomedownloader.py"
python3 - <<'PYEOF'
import sys
sys.path.insert(0, ".")
from harmonizomedownloader import download_datasets
download_datasets(
    [('Achilles Cell Line Gene Essentiality Profiles', 'achilles')],
    ['gene_attribute_edges.txt.gz'],
)
PYEOF
mv -f "Achilles Cell Line Gene Essentiality Profiles/gene_attribute_edges.txt.gz" .
rm -rf "Achilles Cell Line Gene Essentiality Profiles"
echo "[achilles_HMZ] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
