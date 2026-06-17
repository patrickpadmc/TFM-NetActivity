"""
build_entrez_ensg_map.py
Builds Entrez Gene ID -> ENSG mapping via UniProt bridge from Bioteque files.
Saves: data/processed/entrez2ensg.tsv
"""
import csv, sys

GID2UNIPROT  = "/beegfs/home/ppadmoremcc/work/external/bioteque/metadata/mappings/GEN/gid2uniprot.tsv"
ENS2UNIPROT  = "/beegfs/home/ppadmoremcc/work/external/bioteque/metadata/mappings/GEN/ens2uniprot.tsv"
OUT          = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed/entrez2ensg.tsv"

# UniProt -> ENSG
print("Cargando ens2uniprot...")
uniprot2ensg = {}
with open(ENS2UNIPROT) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            ensg, uniprot = parts
            # one UniProt can map to multiple ENSGs; keep first (canonical)
            if uniprot not in uniprot2ensg:
                uniprot2ensg[uniprot] = ensg
print(f"  UniProt->ENSG entries: {len(uniprot2ensg):,}")

# Entrez -> UniProt -> ENSG
print("Cargando gid2uniprot y cruzando...")
entrez2ensg = {}
unmapped = 0
with open(GID2UNIPROT) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            entrez, uniprot = parts
            ensg = uniprot2ensg.get(uniprot)
            if ensg:
                entrez2ensg[entrez] = ensg
            else:
                unmapped += 1

print(f"  Entrez->ENSG mapeados: {len(entrez2ensg):,}")
print(f"  Entrez sin ENSG:       {unmapped:,}")

with open(OUT, "w") as f:
    f.write("entrez_gene_id\tensembl_gene_id\n")
    for entrez, ensg in sorted(entrez2ensg.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 0):
        f.write(f"{entrez}\t{ensg}\n")

print(f"Guardado: {OUT}")
