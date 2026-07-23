#!/bin/bash
#SBATCH --job-name=check_gen_homogenization
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/analysis/logs/%x_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/analysis/logs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

export LC_ALL=C
GRAPH=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/integrated/graph_v2/graph_edges.tsv
TMPSORT=/beegfs/home/ppadmoremcc/work/tmp_sort
mkdir -p "$TMPSORT"

awk -F'\t' 'NR>1 { if ($3=="GEN") print $1; if ($4=="GEN") print $2 }' "$GRAPH" \
    | sort -T "$TMPSORT" --parallel="$SLURM_CPUS_PER_TASK" -u > /tmp/gen_ids_v2.txt

echo "Total GEN unicos:"
wc -l /tmp/gen_ids_v2.txt

echo "GEN que NO empiezan con ENSG:"
grep -cv '^ENSG' /tmp/gen_ids_v2.txt

echo "Muestra de los que no son ENSG (si hay):"
grep -v '^ENSG' /tmp/gen_ids_v2.txt | head -20