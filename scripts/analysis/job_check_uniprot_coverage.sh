#!/bin/bash
#SBATCH --job-name=check_uniprot_coverage
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

python TFM-NetActivity/scripts/analysis/check_uniprot_coverage.py \
    --graph-edges TFM-NetActivity/data/processed/integrated/graph_v2/graph_edges.tsv \
    --ens2uniprot external/bioteque/metadata/mappings/GEN/ens2uniprot.tsv \
    --embeddings-h5 TFM-NetActivity/data/raw/databases/uniprot_embeddings/UP000005640_9606/per-protein.h5 \
    --out-dir TFM-NetActivity/data/processed/analysis/uniprot_coverage