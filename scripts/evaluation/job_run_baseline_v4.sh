#!/bin/bash
#SBATCH --job-name=baseline_v4
#SBATCH --output=scripts/evaluation/logs/baseline_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/baseline_v4_%j.out
#SBATCH --partition=short
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/evaluation/run_baseline_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out data/processed/evaluation/graph_v4/baseline_results_v4.tsv
