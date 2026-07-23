#!/bin/bash
#SBATCH --job-name=spectral_emb
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --mem=220G
#SBATCH --cpus-per-task=48
#SBATCH --output=scripts/embeddings/spectral/logs/spectral_emb_%j.out
#SBATCH --error=scripts/embeddings/spectral/logs/spectral_emb_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -e

cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p data/processed/embeddings/spectral
mkdir -p scripts/embeddings/spectral/logs

python scripts/embeddings/spectral/spectral_embedding.py \
    --edges data/processed/integrated/graph_v2/graph_edges.tsv \
    --out-dir data/processed/embeddings/spectral \
    --k 256