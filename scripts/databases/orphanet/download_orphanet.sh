#!/bin/bash
#SBATCH --job-name=download_orphanet
#SBATCH --output=logs/download_orphanet_%j.out
#SBATCH --error=logs/download_orphanet_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

OUTDIR="/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/orphanet"
mkdir -p "${OUTDIR}"
OUTFILE="${OUTDIR}/en_product6.xml"

echo "[$(date)] Descargando Orphanet en_product6.xml..."

wget --retry-connrefused --waitretry=5 --tries=5 \
     --show-progress \
     -O "${OUTFILE}" \
     "https://www.orphadata.com/data/xml/en_product6.xml"

echo "[$(date)] Verificando que el XML no está vacío y termina correctamente..."
if [[ ! -s "${OUTFILE}" ]]; then
    echo "[$(date)] ERROR: fichero vacío."
    exit 1
fi

LAST_LINE=$(tail -1 "${OUTFILE}")
if [[ "${LAST_LINE}" == *"</JDBOR>"* ]]; then
    echo "[$(date)] XML completo y bien formado."
else
    echo "[$(date)] ADVERTENCIA: el XML puede estar incompleto. Última línea: ${LAST_LINE}"
fi

echo "[$(date)] Tamaño final:"
ls -lh "${OUTFILE}"
