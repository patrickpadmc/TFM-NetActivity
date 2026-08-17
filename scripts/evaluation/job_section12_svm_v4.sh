#!/bin/bash
#SBATCH --job-name=svm_s12
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/logs/svm_s12_%x_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/logs/svm_s12_%x_%j.err

# Uso: sbatch --job-name=svm_<nombre> [--mem=64G] job_svm_generic_v4.sh <source_name> <source_path> [--raw-external]
set -e
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
source .venv/bin/activate

SOURCE_NAME=$1
SOURCE_PATH=$2
EXTRA_ARGS=$3

python3 scripts/evaluation/run_section12_svm.py \
    --source-name "$SOURCE_NAME" \
    --source-path "$SOURCE_PATH" \
    --eval-dir data/processed/evaluation/graph_v4 \
    --out-dir data/processed/analysis/seccion12_svm_link_prediction \
    $EXTRA_ARGS
