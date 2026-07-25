#!/bin/bash
#SBATCH --job-name=opentargets_process
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/opentargets/logs/%x_%j.out
#SBATCH --error=scripts/databases/opentargets/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/opentargets/script.py \
  --raw-dir data/raw/databases/opentargets \
  --out-dir data/processed/databases/opentargets \
  
