#!/bin/bash
#SBATCH --job-name=n2v_node_strict_v4
#SBATCH --output=scripts/evaluation/logs/n2v_node_strict_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/n2v_node_strict_v4_%j.err
#SBATCH --partition=medium
#SBATCH --nodelist=nodo07
#SBATCH --mem=150G
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

mkdir -p data/processed/evaluation/graph_v4/node2vec_final

cd scripts/evaluation
python3 run_node2vec_node_strict_v4.py \
    --train-graph ../../data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --eval-dir ../../data/processed/evaluation/graph_v4 \
    --p 2.0 --q 0.25 \
    --out ../../data/processed/evaluation/graph_v4/node2vec_final/node2vec_results_node_strict.tsv \
    --workers 24
