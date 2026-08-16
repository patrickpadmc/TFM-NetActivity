#!/bin/bash
#SBATCH --job-name=go_f1
#SBATCH --partition=short
#SBATCH --time=09:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/logs/go_f1_%A_%a.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/logs/go_f1_%A_%a.err

set -e
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
source .venv/bin/activate

NAMES=(
  spectral node2vec hin2vec gae_structural gat_structural
  gae_uniprot_fixed gae_uniprot_finetuned
  gat_uniprot_fixed gat_uniprot_finetuned
  uniprot_raw genept_raw
)
FILES=(
  spectral_structural_seed42.tsv node2vec_structural_seed42.tsv hin2vec_structural_seed42.tsv
  gae_structural_seed42.tsv gat_structural_seed42.tsv
  gae_uniprot_fixed_seed42.tsv gae_uniprot_finetuned_seed42.tsv
  gat_uniprot_fixed_seed42.tsv gat_uniprot_finetuned_seed42.tsv
  uniprot_raw_external.tsv genept_raw_external.tsv
)

python3 scripts/evaluation/run_go_benchmark.py \
  --source-name "${NAMES[$SLURM_ARRAY_TASK_ID]}" \
  --source-path "data/processed/analysis/embedding_catalog_v4/${FILES[$SLURM_ARRAY_TASK_ID]}" \
  --data-path data/processed/evaluation/functional_benchmark_v1/go_reference_splits_common_ensg.tsv \
  --out-dir data/processed/analysis/functional_benchmark_v1/go
