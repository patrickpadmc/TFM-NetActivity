#!/bin/bash
#SBATCH --job-name=common_subgraph_v4
#SBATCH --output=scripts/evaluation/logs/common_subgraph_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/common_subgraph_v4_%j.out
#SBATCH --partition=short
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/evaluation/build_common_subgraph_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --out data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
    --report data/processed/evaluation/graph_v4/common_subgraph_report.json
