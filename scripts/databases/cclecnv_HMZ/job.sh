#!/bin/bash
#SBATCH --job-name=cclecnv_HMZ_process
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/cclecnv_HMZ/logs/%x_%j.out
#SBATCH --error=scripts/databases/cclecnv_HMZ/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/cclecnv_HMZ/script.py \
  --raw-dir data/raw/databases/cclecnv_HMZ \
  --out-dir data/processed/databases/cclecnv_HMZ \
  --symbol2ensg data/metadata/symbol2ensg.tsv --model-csv data/raw/databases/cclemut_HMZ/Model.csv
