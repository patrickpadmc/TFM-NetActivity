#!/bin/bash
#SBATCH --job-name=node2vec_v3_convert
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/node2vec/logs/node2vec_v3_convert_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/node2vec/logs/node2vec_v3_convert_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

PROJECT_ROOT=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
OUT_DIR=${PROJECT_ROOT}/data/processed/embeddings/node2vec_v3
EMB_OUT=${OUT_DIR}/node2vec_v3.emb

cd ${PROJECT_ROOT}
source .venv/bin/activate

python scripts/embeddings/node2vec/convert_embeddings_v3.py \
    --emb ${EMB_OUT} \
    --out-dir ${OUT_DIR} \
    --dim 256
