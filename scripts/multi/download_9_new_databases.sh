#!/bin/bash
#SBATCH --job-name=download_9_new_dbs
#SBATCH --output=logs/download_9_new_dbs_%j.out
#SBATCH --error=logs/download_9_new_dbs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

BASE="/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

check_nonempty() {
    local f="$1"
    if [[ ! -s "$f" ]]; then
        log "ERROR: $f está vacío o no existe. Abortando."
        exit 1
    fi
    log "OK: $(du -sh "$f" | cut -f1)  $f"
}

# Descarga solo si el fichero no existe o está vacío
wget_download() {
    local out="$1" url="$2"
    if [[ -s "$out" ]]; then
        log "SKIP (ya existe): $out"
        return
    fi
    wget --retry-connrefused --waitretry=10 --tries=5 \
         --timeout=120 --show-progress -q \
         -O "$out" "$url"
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. HPA v25.0 — fichero unificado (RNA + Proteome)
# ─────────────────────────────────────────────────────────────────────────────
log "=== 1/9  HPA v25.0 ==="
mkdir -p "${BASE}/hpa"
wget_download "${BASE}/hpa/proteinatlas.tsv.zip" \
    "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
if [[ ! -s "${BASE}/hpa/proteinatlas.tsv" ]]; then
    unzip -o "${BASE}/hpa/proteinatlas.tsv.zip" -d "${BASE}/hpa/"
fi
check_nonempty "${BASE}/hpa/proteinatlas.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# 2. Jensen TISSUES curated (actualización semanal)
# ─────────────────────────────────────────────────────────────────────────────
log "=== 2/9  Jensen TISSUES curated ==="
mkdir -p "${BASE}/jensen_tissues"
wget_download "${BASE}/jensen_tissues/human_tissue_knowledge_full.tsv" \
    "https://download.jensenlab.org/human_tissue_knowledge_full.tsv"
check_nonempty "${BASE}/jensen_tissues/human_tissue_knowledge_full.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# 3. STRING v12.0 — Homo sapiens
# ─────────────────────────────────────────────────────────────────────────────
log "=== 3/9  STRING v12.0 ==="
mkdir -p "${BASE}/string"
wget_download "${BASE}/string/9606.protein.links.detailed.v12.0.txt.gz" \
    "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz"
wget_download "${BASE}/string/9606.protein.aliases.v12.0.txt.gz" \
    "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
gzip -t "${BASE}/string/9606.protein.links.detailed.v12.0.txt.gz" && log "OK: STRING links" || { log "ERROR: STRING links corrupto"; exit 1; }
gzip -t "${BASE}/string/9606.protein.aliases.v12.0.txt.gz"        && log "OK: STRING aliases" || { log "ERROR: STRING aliases corrupto"; exit 1; }

# ─────────────────────────────────────────────────────────────────────────────
# 4. OmniPath — red de señalización
# ─────────────────────────────────────────────────────────────────────────────
log "=== 4/9  OmniPath ==="
mkdir -p "${BASE}/omnipath"
wget_download "${BASE}/omnipath/omnipath_interactions.tsv" \
    "https://omnipathdb.org/interactions?datasets=omnipath&format=tsv&genesymbols=1"
check_nonempty "${BASE}/omnipath/omnipath_interactions.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# 5. DoRothEA AB — red regulatoria TF-gen (tiers A+B)
# ─────────────────────────────────────────────────────────────────────────────
log "=== 5/9  DoRothEA AB ==="
mkdir -p "${BASE}/dorothea"
wget_download "${BASE}/dorothea/dorothea_AB_interactions.tsv" \
    "https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B&format=tsv&genesymbols=1"
check_nonempty "${BASE}/dorothea/dorothea_AB_interactions.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# 6. Reactome — mapeo ENSG -> vías
# ─────────────────────────────────────────────────────────────────────────────
log "=== 6/9  Reactome ==="
mkdir -p "${BASE}/reactome"
wget_download "${BASE}/reactome/Ensembl2Reactome_All_Levels.txt" \
    "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
wget_download "${BASE}/reactome/ReactomePathways.txt" \
    "https://reactome.org/download/current/ReactomePathways.txt"
check_nonempty "${BASE}/reactome/Ensembl2Reactome_All_Levels.txt"
check_nonempty "${BASE}/reactome/ReactomePathways.txt"

# ─────────────────────────────────────────────────────────────────────────────
# 7. BioGRID v4.4.244 — Homo sapiens (URL directa con versión)
# ─────────────────────────────────────────────────────────────────────────────
log "=== 7/9  BioGRID v4.4.244 ==="
mkdir -p "${BASE}/biogrid"
wget_download "${BASE}/biogrid/BIOGRID-Homo_sapiens-4.4.244.tab3.zip" \
    "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.244/BIOGRID-ORGANISM-Homo_sapiens-4.4.244.tab3.zip"
if ! ls "${BASE}/biogrid/"BIOGRID-ORGANISM-Homo_sapiens-*.tab3.txt &>/dev/null; then
    unzip -o "${BASE}/biogrid/BIOGRID-Homo_sapiens-4.4.244.tab3.zip" -d "${BASE}/biogrid/"
fi
BIOGRID_FILE=$(ls "${BASE}/biogrid/"BIOGRID-ORGANISM-Homo_sapiens-*.tab3.txt 2>/dev/null | head -1)
[[ -z "$BIOGRID_FILE" ]] && { log "ERROR: fichero BioGRID no encontrado tras descomprimir"; exit 1; }
check_nonempty "$BIOGRID_FILE"

# ─────────────────────────────────────────────────────────────────────────────
# 8. GWAS Catalog — todas las asociaciones (GRCh38)
# ─────────────────────────────────────────────────────────────────────────────
log "=== 8/9  GWAS Catalog ==="
mkdir -p "${BASE}/gwas_catalog"
wget_download "${BASE}/gwas_catalog/gwas_catalog_associations.tsv" \
    "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
check_nonempty "${BASE}/gwas_catalog/gwas_catalog_associations.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# Resumen final
# ─────────────────────────────────────────────────────────────────────────────
log "=== DESCARGA COMPLETADA ==="
du -sh "${BASE}/hpa/"* "${BASE}/jensen_tissues/"* "${BASE}/string/"* \
       "${BASE}/omnipath/"* "${BASE}/dorothea/"* "${BASE}/reactome/"* \
       "${BASE}/biogrid/"* "${BASE}/gwas_catalog/"*

