#!/bin/bash
#SBATCH --job-name=opentargets_dl
#SBATCH --output=scripts/databases/opentargets/logs/%x_%j.out
#SBATCH --error=scripts/databases/opentargets/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
bash scripts/databases/opentargets/get_data.sh "$PWD/data/raw/databases/opentargets"
