#!/bin/bash
#SBATCH --job-name=spectral_subgraph_v4
#SBATCH --output=scripts/evaluation/logs/spectral_subgraph_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/spectral_subgraph_v4_%j.out
#SBATCH --partition=medium
#SBATCH --nodelist=nodo07
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/evaluation/spectral_embedding_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
    --out data/processed/evaluation/graph_v4/embeddings/spectral_v4_subgraph \
    --dim 128
