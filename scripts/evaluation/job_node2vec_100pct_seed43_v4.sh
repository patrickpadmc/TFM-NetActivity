#!/bin/bash
#SBATCH --job-name=n2v_100p_s43_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/n2v_100pct_s43_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/n2v_100pct_s43_v4_%j.err
#SBATCH --partition=medium
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
mkdir -p data/processed/evaluation/graph_v4/benchmark_multiscale/node2vec/100pct_seed43
python3 scripts/evaluation/run_node2vec_final_v4.py \
    --train-graph data/processed/evaluation/graph_v4/benchmark_multiscale/graph_benchmark_100pct_v4.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --best-pq data/processed/evaluation/graph_v4/node2vec_grid/best_pq_per_task.tsv \
    --out-dir data/processed/evaluation/graph_v4/benchmark_multiscale/node2vec/100pct_seed43 \
    --workers 16 \
    --seed 43
