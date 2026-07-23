#!/bin/bash
#SBATCH --job-name=find_isolated_nodes
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/find_isolated_nodes_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/embeddings/spectral/logs/find_isolated_nodes_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -e

cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity

LABELS=data/processed/embeddings/spectral/connectivity_summary_node_labels.tsv
EDGES=data/processed/integrated/graph_v2/graph_edges.tsv
OUT_EDGES=data/processed/embeddings/spectral/isolated_nodes_edges.tsv
OUT_IDS=data/processed/embeddings/spectral/isolated_node_ids.tsv
TMP_IDS=/tmp/isolated_ids_${SLURM_JOB_ID}.txt

echo "Identificando la etiqueta de la componente gigante..."
giant=$(tail -n +2 "$LABELS" | cut -f2 | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
echo "Componente gigante: $giant"

echo "Extrayendo nodos fuera de la componente gigante..."
tail -n +2 "$LABELS" | awk -F'\t' -v g="$giant" '$2 != g {print $1}' > "$TMP_IDS"
cp "$TMP_IDS" "$OUT_IDS"
echo "Nodos aislados encontrados:"
cat "$TMP_IDS"

echo "Buscando estas filas en graph_edges.tsv (109M filas, puede tardar)..."
awk -F'\t' '
NR==FNR { ids[$1]=1; next }
FNR==1 { print; next }
($1 in ids) || ($2 in ids) { print }
' "$TMP_IDS" "$EDGES" > "$OUT_EDGES"

n_matches=$(($(wc -l < "$OUT_EDGES") - 1))
echo "Filas encontradas: $n_matches"
echo "Resultado guardado en: $OUT_EDGES"

rm -f "$TMP_IDS"
echo "Listo"