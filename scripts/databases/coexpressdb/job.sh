#!/bin/bash
#SBATCH --job-name=coexpressdb_process
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/coexpressdb/logs/%x_%j.out
#SBATCH --error=scripts/databases/coexpressdb/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/coexpressdb/script.py \
  --raw-dir data/raw/databases/coexpressdb \
  --out-dir data/processed/databases/coexpressdb \
  --entrez2ensg data/metadata/entrez2ensg.tsv
