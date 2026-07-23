#!/bin/bash
#SBATCH --job-name=count_nodes_edges
#SBATCH --partition=short
#SBATCH --time=00:15:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=scripts/integrated/graph/logs/%x_%j.out
#SBATCH --error=scripts/integrated/graph/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate

cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

python3 scripts/integrated/count_nodes_edges_by_type.py
