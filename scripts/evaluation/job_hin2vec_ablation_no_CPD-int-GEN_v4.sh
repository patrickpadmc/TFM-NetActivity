#!/bin/bash
#SBATCH --job-name=hin2vec_abl_CPD-int-GEN_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/hin2vec_ablation_no_CPD-int-GEN_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/hin2vec_ablation_no_CPD-int-GEN_v4_%j.err
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
mkdir -p data/processed/evaluation/graph_v4/benchmark_multiscale/ablation/no_CPD-int-GEN
python3 scripts/evaluation/run_hin2vec_v4.py \
    --edges data/processed/evaluation/graph_v4/benchmark_multiscale/ablation/graph_ablation_no_CPD-int-GEN_v4.tsv \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out-dir data/processed/evaluation/graph_v4/benchmark_multiscale/ablation/no_CPD-int-GEN \
    --dim 128 \
    --walk-length 6 \
    --max-epochs 200 \
    --patience 20 \
    --check-every 5 \
    --seed 42
