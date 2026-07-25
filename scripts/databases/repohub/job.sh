#!/bin/bash
#SBATCH --job-name=repohub_process
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/repohub/logs/%x_%j.out
#SBATCH --error=scripts/databases/repohub/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/repohub/script.py \
  --raw-dir data/raw/databases/repohub \
  --out-dir data/processed/databases/repohub \
  --symbol2ensg data/metadata/symbol2ensg.tsv
