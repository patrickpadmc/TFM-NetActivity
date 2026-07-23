#!/bin/bash
#SBATCH --job-name=ctdchemgen_graph_v2
#SBATCH --partition=short
#SBATCH --time=00:45:00
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

echo "================================================"
echo "[$(date)] PASO 1 -- Procesar ctdchemgen (3 relaciones expression)"
echo "================================================"
python3 scripts/databases/ctdchemgen/script.py

echo ""
echo "================================================"
echo "[$(date)] PASO 2 -- Construir grafo v2 (17 + ctdchemgen)"
echo "================================================"
python3 scripts/integrated/build_graph_table_v2.py

echo ""
echo "================================================"
echo "[$(date)] COMPLETADO."
echo "================================================"
