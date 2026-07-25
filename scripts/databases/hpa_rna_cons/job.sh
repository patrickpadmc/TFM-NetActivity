#!/bin/bash
#SBATCH --job-name=hpa_rna_cons_process
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/hpa_rna_cons/logs/%x_%j.out
#SBATCH --error=scripts/databases/hpa_rna_cons/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/hpa_rna_cons/script.py \
  --raw-dir data/raw/databases/hpa_rna_cons \
  --out-dir data/processed/databases/hpa_rna_cons \
  
