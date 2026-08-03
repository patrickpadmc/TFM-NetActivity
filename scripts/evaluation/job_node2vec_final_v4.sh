#!/bin/bash
#SBATCH --job-name=n2v_final_v4
#SBATCH --output=scripts/evaluation/logs/n2v_final_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/n2v_final_v4_%j.err
#SBATCH --partition=medium
#SBATCH --nodelist=nodo11
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

mkdir -p data/processed/evaluation/graph_v4/node2vec_final

python3 scripts/evaluation/run_node2vec_final_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --best-pq data/processed/evaluation/graph_v4/node2vec_grid/best_pq_per_task.tsv \
    --out-dir data/processed/evaluation/graph_v4/node2vec_final \
    --workers 32
