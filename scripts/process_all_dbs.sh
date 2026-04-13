#!/bin/bash
#SBATCH --job-name=process_all_dbs
#SBATCH --output=logs/process_all_dbs_%j.out
#SBATCH --error=logs/process_all_dbs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate

SCRIPTS=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts
PROCESSED=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed
mkdir -p "${PROCESSED}"

echo "================================================"
echo "[$(date)] PASO 1 — Mapeo Entrez -> ENSG"
echo "================================================"
python3 "${SCRIPTS}/build_entrez_ensg_map.py"

echo ""
echo "================================================"
echo "[$(date)] PASO 2 — Procesar CTD"
echo "================================================"
python3 "${SCRIPTS}/process_ctd.py"

echo ""
echo "================================================"
echo "[$(date)] PASO 3 — Procesar HPO"
echo "================================================"
python3 "${SCRIPTS}/process_hpo.py"

echo ""
echo "================================================"
echo "[$(date)] PASO 4 — Procesar Orphanet"
echo "================================================"
python3 "${SCRIPTS}/process_orphanet.py"

echo ""
echo "================================================"
echo "[$(date)] PASO 5 — Dashboard Open Targets v2"
echo "================================================"
python3 "${SCRIPTS}/make_dashboard_ot_26.03_v2.py"

echo ""
echo "================================================"
echo "[$(date)] COMPLETADO. Ficheros generados:"
echo "================================================"
ls -lh "${PROCESSED}/"
