#!/bin/bash
#SBATCH --job-name=gae_v4
#SBATCH --output=scripts/evaluation/logs/gae_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/gae_v4_%j.out
#SBATCH --partition=short
#SBATCH --gres=gpu:tesla_p100:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

mkdir -p data/processed/evaluation/graph_v4/gae

python3 scripts/evaluation/run_gae_v4.py \
    --train-graph data/processed/evaluation/graph_v4/graph_edges_train_only.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out-dir data/processed/evaluation/graph_v4/gae \
    --dim 128 \
    --max-epochs 200 \
    --patience 10
