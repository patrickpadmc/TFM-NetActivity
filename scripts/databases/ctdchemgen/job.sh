#!/bin/bash
#SBATCH --job-name=ctdchemgen
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=scripts/databases/ctdchemgen/logs/%x_%j.out
#SBATCH --error=scripts/databases/ctdchemgen/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

cd ~/work/TFM-NetActivity

python3 scripts/databases/ctdchemgen/script.py \
  --raw-dir data/raw/databases/ctdchemgen \
  --out-dir data/processed/databases/ctdchemgen \
  --entrez-mapping data/metadata/entrez2ensg.tsv