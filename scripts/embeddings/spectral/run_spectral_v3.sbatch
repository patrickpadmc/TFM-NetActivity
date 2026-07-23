#!/bin/bash
#SBATCH --job-name=spectral_v3
#SBATCH --partition=short
#SBATCH --mem=120G
#SBATCH --cpus-per-task=48
#SBATCH --time=04:00:00
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/spectral_v3_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/spectral_v3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

PROJECT_ROOT=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
EDGES=${PROJECT_ROOT}/data/processed/integrated/graph_v3/graph_edges.tsv
OUT_DIR=${PROJECT_ROOT}/data/processed/embeddings/spectral_v3

cd ${PROJECT_ROOT}
source .venv/bin/activate

python scripts/embeddings/spectral/spectral_embedding_v3.py \
    --edges ${EDGES} \
    --out-dir ${OUT_DIR} \
    --k 256
