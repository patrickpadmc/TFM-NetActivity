#!/bin/bash
#SBATCH --job-name=check_connectivity
#SBATCH --partition=medium
#SBATCH --time=02:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=16
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/check_connectivity_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/check_connectivity_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

set -e

cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate

python scripts/embeddings/spectral/check_connectivity.py \
    --edges data/processed/integrated/graph_v2/graph_edges.tsv \
    --out data/processed/embeddings/spectral/connectivity_summary.tsv