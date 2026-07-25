#!/bin/bash
#SBATCH --job-name=conn_comp_v4
#SBATCH --output=scripts/analysis/logs/%x_%j.out
#SBATCH --error=scripts/analysis/logs/%x_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/analysis/connected_components_v4.py \
    --edges data/processed/integrated/graph_v4/graph_edges.tsv \
    --out-dir data/processed/integrated/graph_v4/connectivity
