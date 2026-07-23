#!/bin/bash
#SBATCH --job-name=graph_v2_stats_flag_mappings
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=TFM-NetActivity/scripts/analysis/logs/%x_%j.out
#SBATCH --error=TFM-NetActivity/scripts/analysis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate

cd /beegfs/home/ppadmoremcc/work

python TFM-NetActivity/scripts/analysis/graph_v2_stats.py \
    --graph-edges TFM-NetActivity/data/processed/integrated/graph_v2/graph_edges.tsv \
    --out-dir TFM-NetActivity/data/processed/analysis/graph_v2_stats

python TFM-NetActivity/scripts/analysis/flag_ambiguous_mappings.py \
    --mapping-full TFM-NetActivity/data/processed/analysis/uniprot_coverage/ensg_uniprot_mapping_full.tsv \
    --out-dir TFM-NetActivity/data/processed/analysis/uniprot_coverage