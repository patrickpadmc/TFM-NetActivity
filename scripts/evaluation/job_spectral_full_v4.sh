#!/bin/bash
#SBATCH --job-name=spectral_full_v4
#SBATCH --output=scripts/evaluation/logs/spectral_full_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/spectral_full_v4_%j.out
#SBATCH --partition=medium
#SBATCH --nodelist=nodo07
#SBATCH --mem=450G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

mkdir -p data/processed/evaluation/graph_v4/embeddings

python3 scripts/evaluation/spectral_embedding_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --out data/processed/evaluation/graph_v4/embeddings/spectral_v4_full \
    --dim 128
