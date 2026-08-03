#!/bin/bash
#SBATCH --job-name=spectral_eval_v4
#SBATCH --output=scripts/evaluation/logs/spectral_eval_v4_%j.out
#SBATCH --error=scripts/evaluation/logs/spectral_eval_v4_%j.err
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/evaluation/run_spectral_eval_v4.py \
    --embeddings data/processed/evaluation/graph_v4/embeddings/spectral_v4_subgraph.npz \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out data/processed/evaluation/graph_v4/spectral_results.tsv
