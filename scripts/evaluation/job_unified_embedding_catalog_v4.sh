#!/bin/bash
#SBATCH --job-name=embedding_catalog_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/embedding_catalog_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/evaluation/logs/embedding_catalog_v4_%j.err
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
python3 scripts/evaluation/build_unified_embedding_catalog_v4.py \
    --subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
    --out-dir data/processed/analysis/embedding_catalog_v4
