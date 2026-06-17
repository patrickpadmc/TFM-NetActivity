#!/bin/bash
#SBATCH --job-name=download_ctd
#SBATCH --output=logs/download_ctd_%j.out
#SBATCH --error=logs/download_ctd_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

OUTDIR="/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/ctd/"
BASE_URL="https://ctdbase.org/reports"

echo "[$(date)] Iniciando descarga de CTD..."

# Asociaciones gen-enfermedad (curadas + inferidas)
wget -q --show-progress -O "${OUTDIR}/CTD_genes_diseases.tsv.gz" \
    "${BASE_URL}/CTD_genes_diseases.tsv.gz"

echo "[$(date)] Descarga completada."
echo "Ficheros descargados:"
ls -lh "${OUTDIR}/"
