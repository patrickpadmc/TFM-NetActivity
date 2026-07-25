#!/bin/bash
#SBATCH --job-name=ctdchemis_process
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=scripts/databases/ctdchemis/logs/%x_%j.out
#SBATCH --error=scripts/databases/ctdchemis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd ~/work/TFM-NetActivity

python3 scripts/databases/ctdchemis/script.py \
  --raw-dir data/raw/databases/ctdchemis \
  --out-dir data/processed/databases/ctdchemis \
  --ctd-mapping /home/ppadmoremcc/work/external/bioteque/metadata/mappings/CPD/ctd.tsv
