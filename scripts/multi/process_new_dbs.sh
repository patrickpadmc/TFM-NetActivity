#!/bin/bash
#SBATCH --job-name=process_new_dbs
#SBATCH --output=logs/process_new_dbs_%j.out
#SBATCH --error=logs/process_new_dbs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate

SCRIPTS=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts
mkdir -p /beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed

echo "================================================"
echo "[$(date)] 1/6 — HPA"
echo "================================================"
python3 "${SCRIPTS}/process_hpa.py"

echo ""
echo "================================================"
echo "[$(date)] 2/6 — Jensen TISSUES"
echo "================================================"
python3 "${SCRIPTS}/process_jensen.py"

echo ""
echo "================================================"
echo "[$(date)] 3/6 — STRING"
echo "================================================"
python3 "${SCRIPTS}/process_string.py"

echo ""
echo "================================================"
echo "[$(date)] 4-5-6/6 — OmniPath + DoRothEA + Reactome"
echo "================================================"
python3 "${SCRIPTS}/process_omnipath_dorothea_reactome.py"

echo ""
echo "================================================"
echo "[$(date)] COMPLETADO — ficheros generados:"
echo "================================================"
ls -lh /beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/*.tsv \
        /beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/dashboard_*.html
