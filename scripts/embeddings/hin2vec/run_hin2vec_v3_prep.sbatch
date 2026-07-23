#!/bin/bash
#SBATCH --job-name=hin2vec_v3_prep
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/hin2vec/logs/hin2vec_v3_prep_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/hin2vec/logs/hin2vec_v3_prep_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

# Correr UNA SOLA VEZ antes de lanzar las 4 variantes (run_hin2vec_v3_*.sbatch).
# Las 4 variantes leen el parquet cacheado que este job genera.

set -euo pipefail

PROJECT_ROOT=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
EDGES=${PROJECT_ROOT}/data/processed/integrated/graph_v3/graph_edges.tsv
OUT_DIR=${PROJECT_ROOT}/data/processed/embeddings/hin2vec_v3
CACHED_PARQUET=${OUT_DIR}/graph_v3_with_edge_class.parquet

cd ${PROJECT_ROOT}
source .venv/bin/activate
mkdir -p ${OUT_DIR}

python scripts/embeddings/hin2vec/prepare_edge_class_v3.py \
    --edges ${EDGES} \
    --out ${CACHED_PARQUET}
