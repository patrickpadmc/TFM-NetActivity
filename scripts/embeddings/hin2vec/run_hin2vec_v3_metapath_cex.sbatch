#!/bin/bash
#SBATCH --job-name=hin2vec_v3_metapath_cex
#SBATCH --partition=long
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:tesla_p100:1
#SBATCH --time=48:00:00
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/hin2vec/logs/hin2vec_v3_metapath_cex_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/hin2vec/logs/hin2vec_v3_metapath_cex_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

# Requiere haber corrido run_hin2vec_v3_prep.sbatch antes.

set -euo pipefail

PROJECT_ROOT=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
OUT_DIR=${PROJECT_ROOT}/data/processed/embeddings/hin2vec_v3
CACHED_PARQUET=${OUT_DIR}/graph_v3_with_edge_class.parquet
METAPATH_RELATIONS=${PROJECT_ROOT}/scripts/embeddings/hin2vec/metapath_relations_v3.tsv
RUN_NAME=metapath_cex
EDGES_PARQUET=${OUT_DIR}/hin_input_v3_${RUN_NAME}.parquet

cd ${PROJECT_ROOT}
source .venv/bin/activate

if [ ! -f "${EDGES_PARQUET}" ]; then
    python scripts/embeddings/hin2vec/build_hin_input_v3.py \
        --cached-parquet ${CACHED_PARQUET} \
        --out ${EDGES_PARQUET} \
        --mode metapath \
        --metapath-relations ${METAPATH_RELATIONS} \
        --include-cex
fi

python scripts/embeddings/hin2vec/train_hin2vec_v3.py \
    --edges-parquet ${EDGES_PARQUET} \
    --out-dir ${OUT_DIR} \
    --run-name ${RUN_NAME} \
    --window 4 \
    --walk 10 \
    --walk-length 100 \
    --embed-size 256 \
    --neg 5 \
    --epochs 10 \
    --batch-size 4096
