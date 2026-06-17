#!/bin/bash
#SBATCH --job-name=tfm_smoke
#SBATCH --output=logs/tfm_smoke_%j.out
#SBATCH --error=logs/tfm_smoke_%j.err
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=short

set -euo pipefail
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
python src/utils/smoke_test.py
