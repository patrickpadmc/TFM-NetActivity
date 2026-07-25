#!/usr/bin/env bash
# =============================================================================
# ctdchemis — get_data.sh
# Fuente: Comparative Toxigenomics Database (CTD)
# Archivos: CTD_chemicals_diseases.tsv.gz
# Contiene: CPD-cau-DIS (compound causes disease) y CPD-trt-DIS (treats disease)
#
# CTD bloquea descargas automatizadas con un captcha de verificacion.
# Este script valida la cabecera gzip del resultado: si CTD devuelve la
# pagina HTML del captcha en vez del fichero real, lo detecta y avisa.
# =============================================================================
set -euo pipefail
OUT_DIR="${1:-$(dirname "$0")/../../data/raw/databases/ctdchemis}"
OUT="$OUT_DIR/CTD_chemicals_diseases.tsv.gz"
URL="https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz"
mkdir -p "$OUT_DIR"

if [ -s "$OUT" ]; then
    echo "Ya existe y no esta vacio: $OUT (omitiendo descarga)"
else
    echo "[ctdchemis] Descargando a $OUT ..."
    wget -O "$OUT" "$URL"
fi

MAGIC=$(head -c 2 "$OUT" | od -An -tx1 | tr -d ' ')
if [ "$MAGIC" != "1f8b" ]; then
    echo "ERROR: el fichero descargado no es un gzip valido."
    echo "CTD probablemente devolvio una pagina de verificacion (captcha)."
    echo "Primeras lineas del fichero:"
    head -c 300 "$OUT"
    echo ""
    echo ""
    echo "Borra el fichero corrupto (rm $OUT) y descargalo manualmente:"
    echo "  1. Abre $URL en el navegador (completa el captcha si aparece)"
    echo "  2. scp el fichero descargado al cluster:"
    echo "     scp CTD_chemicals_diseases.tsv.gz <usuario>@hpclogin.unav.es:$OUT_DIR/"
    exit 1
fi
echo "OK: gzip valido en $OUT"
gzip -t "$OUT" && echo "OK: integridad del gzip verificada"
