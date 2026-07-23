#!/bin/bash
#SBATCH --job-name=build_graph_v3
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/integrated/graph_v3/logs/build_graph_v3_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/integrated/graph_v3/logs/build_graph_v3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -e
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

module load Python/3.11.3-GCCcore-12.3.0
source .venv/bin/activate

python scripts/integrated/graph_v3/build_graph_table_v3.py
