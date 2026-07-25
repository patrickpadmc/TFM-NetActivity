#!/bin/bash
#SBATCH --job-name=download_all_v4
#SBATCH --output=scripts/multi/logs/%x_%j.out
#SBATCH --error=scripts/multi/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -uo pipefail  # sin -e: un fallo en una base no debe abortar el resto

BASE=/beegfs/home/ppadmoremcc/work/TFM-NetActivity
cd "$BASE"
mkdir -p logs/download_v4

FAILED=()

download_one () {
    local db="$1"
    echo "================================================"
    echo "[$(date)] Descargando: $db"
    echo "================================================"
    if bash "scripts/databases/$db/get_data.sh" "$BASE/data/raw/databases/$db" \
        > "logs/download_v4/${db}.log" 2>&1; then
        echo "[$db] OK"
    else
        echo "[$db] FALLO -- ver logs/download_v4/${db}.log"
        FAILED+=("$db")
    fi
}

# omnipath primero: dorothea_AB y dorothea_CD dependen de sus archivos (symlink)
download_one omnipath
download_one dorothea_AB
download_one dorothea_CD

for db in achilles_HMZ cclecnv_HMZ cclemut_HMZ cclerna_HMZ coexpressdb \
          ctdchemgen ctdchemis hpa_proteome hpa_rna_cons jensentissuecurated \
          opentargets reactome repohub string; do
    download_one "$db"
done

echo ""
echo "================================================"
echo "Validando integridad de .gz / .zip descargados"
echo "================================================"
BAD=()
while IFS= read -r -d '' f; do
    case "$f" in
        *.gz)
            gzip -t "$f" 2>/dev/null && echo "OK gzip: $f" || { echo "FALLO gzip (posible captcha/HTML): $f"; BAD+=("$f"); }
            ;;
        *.zip)
            unzip -tq "$f" >/dev/null 2>&1 && echo "OK zip: $f" || { echo "FALLO zip: $f"; BAD+=("$f"); }
            ;;
    esac
done < <(find data/raw/databases -type f \( -name "*.gz" -o -name "*.zip" \) -print0)

echo ""
echo "================================================"
echo "RESUMEN"
echo "================================================"
if [ ${#FAILED[@]} -eq 0 ]; then
    echo "Todas las descargas terminaron sin error de script."
else
    echo "Bases con FALLO en get_data.sh: ${FAILED[*]}"
fi
if [ ${#BAD[@]} -eq 0 ]; then
    echo "Todos los .gz/.zip pasaron la validacion de integridad."
else
    echo "Archivos con integridad SOSPECHOSA (revisar manualmente, posible captcha): ${BAD[*]}"
fi
