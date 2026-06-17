#!/bin/bash
#SBATCH --job-name=parse_orphanet
#SBATCH --output=logs/parse_orphanet_%j.out
#SBATCH --error=logs/parse_orphanet_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es

set -euo pipefail

module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate

python3 - << 'PYEOF'
import xml.etree.ElementTree as ET
import csv

xml_path = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/orphanet/en_product6.xml"
out_path = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/raw/orphanet/orphanet_gene_disease.tsv"

print(f"Parseando {xml_path}...")
tree = ET.parse(xml_path)
root = tree.getroot()

rows = []

for disorder in root.iter("Disorder"):
    orpha_code = disorder.findtext("OrphaCode", default="")
    disease_name = disorder.findtext("Name", default="")

    gene_list = disorder.find("DisorderGeneAssociationList")
    if gene_list is None:
        continue

    for assoc in gene_list.findall("DisorderGeneAssociation"):
        gene = assoc.find("Gene")
        if gene is None:
            continue

        gene_symbol = gene.findtext("Symbol", default="")
        gene_name   = gene.findtext("Name", default="")

        # Entrez Gene ID
        entrez_id = ""
        ensembl_id = ""
        for ext_ref in gene.iter("ExternalReference"):
            source = ext_ref.findtext("Source", default="")
            ref    = ext_ref.findtext("Reference", default="")
            if source == "HGNC":
                pass
            elif source == "Ensembl":
                ensembl_id = ref
            elif source == "OMIM":
                pass

        assoc_type = assoc.findtext(
            "DisorderGeneAssociationType/Name", default="")
        assoc_status = assoc.findtext(
            "DisorderGeneAssociationStatus/Name", default="")

        rows.append({
            "orpha_code":    orpha_code,
            "disease_name":  disease_name,
            "gene_symbol":   gene_symbol,
            "gene_name":     gene_name,
            "ensembl_id":    ensembl_id,
            "association_type":   assoc_type,
            "association_status": assoc_status,
        })

fieldnames = [
    "orpha_code", "disease_name", "gene_symbol",
    "gene_name", "ensembl_id",
    "association_type", "association_status"
]

with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"TSV escrito en: {out_path}")
print(f"Total de asociaciones: {len(rows):,}")
PYEOF
