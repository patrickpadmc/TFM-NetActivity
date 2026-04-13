#!/bin/bash
#SBATCH --job-name=filter_ctd
#SBATCH --output=logs/filter_ctd_%j.out
#SBATCH --error=logs/filter_ctd_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

INDIR="/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/ctd"
INFILE="${INDIR}/CTD_genes_diseases.tsv.gz"
OUTFILE="${INDIR}/CTD_genes_diseases_human.tsv"

echo "[$(date)] Filtrando por Homo sapiens (OrganismID=9606)..."

# Extraer cabecera (líneas que empiezan por '#') + filas con OrganismID 9606
# La columna OrganismID es la columna 9
zcat "${INFILE}" | awk '
    /^#/ { print; next }
    NR == 1 { print; next }
    $9 == "9606" { print }
' > "${OUTFILE}"

echo "[$(date)] Filtrado completado."
echo "Líneas totales en el fichero humano:"
wc -l "${OUTFILE}"
