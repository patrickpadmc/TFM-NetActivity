#!/bin/bash
#SBATCH --job-name=n2v_grid_v4
#SBATCH --output=scripts/evaluation/logs/n2v_grid_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/n2v_grid_v4_%j.out
#SBATCH --partition=medium
#SBATCH --nodelist=nodo11
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/evaluation/run_node2vec_grid_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out-dir data/processed/evaluation/graph_v4/node2vec_grid \
    --workers 32
