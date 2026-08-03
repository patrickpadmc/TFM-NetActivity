#!/bin/bash
#SBATCH --job-name=gae_75pct_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/gae_75pct_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/gae_75pct_v4_%j.err
#SBATCH --partition=short
#SBATCH --gres=gpu:tesla_p100:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
mkdir -p data/processed/evaluation/graph_v4/benchmark_multiscale/gae/75pct
python3 scripts/evaluation/run_gae_v4.py \
    --train-graph data/processed/evaluation/graph_v4/benchmark_multiscale/graph_benchmark_75pct_v4.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out-dir data/processed/evaluation/graph_v4/benchmark_multiscale/gae/75pct \
    --dim 128 \
    --max-epochs 200 \
    --patience 10
