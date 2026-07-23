#!/usr/bin/env bash
# =============================================================================
# cclemut_HMZ — get_data.sh
# Fuente: Harmonizome / CCLE Cell Line Gene Mutation Profiles
# Contiene: CLL-mut-GEN (cell mutation gene)
# =============================================================================
set -euo pipefail

OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/cclemut_HMZ}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

API_SCRIPT="${OUT_DIR}/harmonizomeapi.py"

if [[ ! -f "$API_SCRIPT" ]]; then
    echo "[cclemut_HMZ] Descargando Harmonizome API..."
    wget --continue --show-progress \
        -O "$API_SCRIPT" \
        "https://maayanlab.cloud/Harmonizome/static/harmonizomeapi.py"
fi

echo "[cclemut_HMZ] Obteniendo ruta de descarga via API..."
python3 - <<'PYEOF'
import sys
sys.path.insert(0, ".")
import harmonizomeapi as api
import urllib.request, os

dataset = api.get("dataset", "CCLE+Cell+Line+Gene+Mutation+Profiles")
download_url = dataset.get("downloadUrl") or dataset.get("download_url")

if not download_url:
    print("ERROR: No se encontró downloadUrl en la respuesta de la API.")
    print("Respuesta:", dataset)
    sys.exit(1)

if not download_url.startswith("http"):
    download_url = "https://maayanlab.cloud" + download_url

print(f"[cclemut_HMZ] URL de descarga: {download_url}")
filename = os.path.basename(download_url)
urllib.request.urlretrieve(download_url, filename)
print(f"[cclemut_HMZ] Descargado: {filename}")
PYEOF

echo "[cclemut_HMZ] Listo. Archivos en: $OUT_DIR"
ls -lh "$OUT_DIR"
