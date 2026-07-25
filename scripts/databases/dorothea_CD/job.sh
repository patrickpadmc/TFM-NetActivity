#!/bin/bash
#SBATCH --job-name=dorothea_CD_process
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/dorothea_CD/logs/%x_%j.out
#SBATCH --error=scripts/databases/dorothea_CD/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/dorothea_CD/script.py \
  --raw-dir data/raw/databases/dorothea_CD \
  --out-dir data/processed/databases/dorothea_CD \
  --ensp2ensg data/metadata/ensp2ensg.tsv.gz
