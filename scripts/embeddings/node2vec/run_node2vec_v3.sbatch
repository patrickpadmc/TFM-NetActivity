#!/bin/bash
#SBATCH --job-name=node2vec_v3
#SBATCH --partition=medium
#SBATCH --mem=150G
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/node2vec/logs/node2vec_v3_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/node2vec/logs/node2vec_v3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

PROJECT_ROOT=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
EDGES=${PROJECT_ROOT}/data/processed/integrated/graph_v3/graph_edges.tsv
OUT_DIR=${PROJECT_ROOT}/data/processed/embeddings/node2vec_v3
EDGELIST=${OUT_DIR}/graph_v3_dedup.edg
EMB_OUT=${OUT_DIR}/node2vec_v3.emb

DIM=256
WALK_LENGTH=80
NUM_WALKS=10
WINDOW_SIZE=10
EPOCHS=5
P=1
Q=1

cd ${PROJECT_ROOT}
source .venv/bin/activate
mkdir -p ${OUT_DIR}

# Etapa 1: construir edgelist deduplicado (skip si ya existe)
if [ ! -f "${EDGELIST}" ]; then
    python scripts/embeddings/node2vec/build_edgelist_v3.py \
        --edges ${EDGES} \
        --out ${EDGELIST}
fi

# Etapa 2: node2vec via pecanpy, modo SparseOTF (>10k nodos, densidad ~1.9%)
pecanpy --input ${EDGELIST} \
    --output ${EMB_OUT} \
    --mode SparseOTF \
    --dimensions ${DIM} \
    --walk-length ${WALK_LENGTH} \
    --num-walks ${NUM_WALKS} \
    --window-size ${WINDOW_SIZE} \
    --epochs ${EPOCHS} \
    --p ${P} \
    --q ${Q} \
    --workers 48

# Etapa 3: convertir a .npy + node_order.tsv (misma convencion que spectral)
python scripts/embeddings/node2vec/convert_embeddings_v3.py \
    --emb ${EMB_OUT} \
    --out-dir ${OUT_DIR} \
    --dim ${DIM}
