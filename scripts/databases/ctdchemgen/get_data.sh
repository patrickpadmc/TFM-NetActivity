#!/bin/bash
# get_data.sh - descarga CTD_chem_gene_ixns.tsv.gz
#
# CTD bloquea descargas automatizadas con un captcha de verificacion
# (mismo problema que ya tuvimos con Reactome). Este script intenta
# wget, pero valida la cabecera gzip del resultado: si CTD devuelve la
# pagina HTML del captcha en vez del fichero real, lo detecta y avisa.
#
# Uso: ./get_data.sh [raw_dir]

set -e

RAW_DIR="${1:-$HOME/data/raw/databases/ctdchemgen}"
URL="https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz"
OUT="$RAW_DIR/CTD_chem_gene_ixns.tsv.gz"

mkdir -p "$RAW_DIR"

if [ -s "$OUT" ]; then
    echo "Ya existe y no esta vacio: $OUT (omitiendo descarga)"
else
    echo "Descargando a $OUT ..."
    wget -O "$OUT" "$URL"
fi

# Validar que es un gzip real y no una pagina de captcha
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
    echo "     scp CTD_chem_gene_ixns.tsv.gz <usuario>@hpclogin.unav.es:$RAW_DIR/"
    exit 1
fi

echo "OK: gzip valido en $OUT"
gzip -t "$OUT" && echo "OK: integridad del gzip verificada"
