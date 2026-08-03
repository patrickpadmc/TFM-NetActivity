#!/bin/bash
#SBATCH --job-name=ext_embeddings_v4
#SBATCH --output=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/analysis/logs/ext_embeddings_v4_%j.out
#SBATCH --error=/beegfs/home/ppadmoremcc/work/TFM-NetActivity/scripts/analysis/logs/ext_embeddings_v4_%j.err
#SBATCH --partition=short
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppadmoremcc@alumni.unav.es
module load Python/3.11.3-GCCcore-12.3.0
source /beegfs/home/ppadmoremcc/work/TFM-NetActivity/.venv/bin/activate
cd /beegfs/home/ppadmoremcc/work/TFM-NetActivity
python3 scripts/analysis/build_external_embeddings_v4.py \
    --subgraph data/processed/evaluation/graph_v4/graph_edges_common_subgraph_v4.tsv \
    --ens2uniprot ../external/bioteque/metadata/mappings/GEN/ens2uniprot.tsv \
    --reviewed ../external/bioteque/metadata/mappings/GEN/human_reviewed.tsv \
    --gname2uniprot ../external/bioteque/metadata/mappings/GEN/gname2uniprot.tsv \
    --uniprot-h5 data/raw/databases/uniprot_embeddings/UP000005640_9606/per-protein.h5 \
    --genept-pickle data/raw/databases/genept/GenePT_gene_protein_embedding_model_3_text.pickle \
    --out-dir data/processed/analysis/external_embeddings_v4
