#!/bin/bash
#SBATCH --job-name=spectral_25pct_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/spectral_25pct_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/spectral_25pct_v4_%j.err
#SBATCH --partition=medium
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
python3 scripts/evaluation/spectral_embedding_v4.py \
    --train-graph data/processed/evaluation/graph_v4/benchmark_multiscale/graph_benchmark_25pct_v4.tsv \
    --out data/processed/evaluation/graph_v4/benchmark_multiscale/embeddings/spectral_v4_25pct \
    --dim 128
python3 scripts/evaluation/run_spectral_eval_v4.py \
    --embeddings data/processed/evaluation/graph_v4/benchmark_multiscale/embeddings/spectral_v4_25pct.npz \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out data/processed/evaluation/graph_v4/benchmark_multiscale/spectral/spectral_results_25pct.tsv
