#!/bin/bash

#SBATCH --job-name=download_hpo
#SBATCH --output=logs/download_hpo_%j.out
#SBATCH --error=logs/download_hpo_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

OUTDIR="/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/hpo"
BASE_URL="https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"

echo "[$(date)] Iniciando descarga de HPO..."

# Anotaciones gen-enfermedad (fichero principal para este proyecto)
wget -q --show-progress -O "${OUTDIR}/genes_to_disease.txt" \
    "${BASE_URL}/genes_to_disease.txt"

# Anotaciones fenotípicas completas (enfermedad -> terminos HP:)
wget -q --show-progress -O "${OUTDIR}/phenotype.hpoa" \
    "${BASE_URL}/phenotype.hpoa"

# Ontología completa en formato OBO (opcional, útil para análisis semántico)
wget -q --show-progress -O "${OUTDIR}/hp.obo" \
    "${BASE_URL}/hp.obo"

echo "[$(date)] Descarga completada."
echo "Ficheros descargados:"
ls -lh "${OUTDIR}/"
