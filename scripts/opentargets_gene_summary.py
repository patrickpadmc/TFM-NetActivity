#!/usr/bin/env python3
import os
import json
import csv
from collections import defaultdict
from statistics import mean

# -----------------------------
# Paths
# -----------------------------
INPUT_DIR = "/beegfs/home/ppadmoremcc/work/external/bioteque/datasets/opentargets/associationByOverallIndirect"
OUT_DIR = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/processed"
META_DIR = "/beegfs/home/ppadmoremcc/work/TFM-NetActivity/data/metadata"

SUMMARY_TSV = os.path.join(OUT_DIR, "opentargets_gene_summary.tsv")
GENERAL_TXT = os.path.join(META_DIR, "opentargets_general_summary.txt")
METADATA_TXT = os.path.join(META_DIR, "opentargets_metadata.txt")

# -----------------------------
# Storage
# -----------------------------
gene_stats = defaultdict(lambda: {
    "n_associations": 0,
    "scores": [],
    "evidence_counts": []
})

total_rows = 0
bad_lines = 0

# -----------------------------
# Read all json files
# -----------------------------
files = sorted([
    f for f in os.listdir(INPUT_DIR)
    if f.endswith(".json")
])

for fname in files:
    fpath = os.path.join(INPUT_DIR, fname)
    with open(fpath, "r") as fh:
        for line in fh:
            try:
                x = json.loads(line)
                gene_id = x["targetId"]
                score = float(x["score"])
                evidence_count = int(x["evidenceCount"])

                gene_stats[gene_id]["n_associations"] += 1
                gene_stats[gene_id]["scores"].append(score)
                gene_stats[gene_id]["evidence_counts"].append(evidence_count)

                total_rows += 1
            except Exception:
                bad_lines += 1

# -----------------------------
# Write per-gene summary
# -----------------------------
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(META_DIR, exist_ok=True)

with open(SUMMARY_TSV, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow([
        "targetId",
        "n_associations",
        "mean_score",
        "max_score",
        "sum_evidenceCount",
        "mean_evidenceCount"
    ])

    for gene_id in sorted(gene_stats.keys()):
        scores = gene_stats[gene_id]["scores"]
        evids = gene_stats[gene_id]["evidence_counts"]

        writer.writerow([
            gene_id,
            gene_stats[gene_id]["n_associations"],
            round(mean(scores), 6),
            round(max(scores), 6),
            sum(evids),
            round(mean(evids), 6)
        ])

# -----------------------------
# Write general summary
# -----------------------------
with open(GENERAL_TXT, "w") as out_f:
    out_f.write("dataset\tOpen Targets\n")
    out_f.write("release\t22.04\n")
    out_f.write("source_path\t%s\n" % INPUT_DIR)
    out_f.write("total_associations\t%d\n" % total_rows)
    out_f.write("unique_genes\t%d\n" % len(gene_stats))
    out_f.write("bad_lines\t%d\n" % bad_lines)

# -----------------------------
# Write metadata
# -----------------------------
with open(METADATA_TXT, "w") as out_f:
    out_f.write("dataset_name\topentargets\n")
    out_f.write("source\tOpen Targets\n")
    out_f.write("release_used\t22.04\n")
    out_f.write("download_source\tftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.04/output/etl/json/\n")
    out_f.write("download_date\t2026-03-15\n")
    out_f.write("gene_id_field\ttargetId\n")
    out_f.write("gene_id_type\tEnsembl Gene ID\n")
    out_f.write("disease_id_field\tdiseaseId\n")
    out_f.write("score_field\tscore\n")
    out_f.write("evidence_count_field\tevidenceCount\n")
    out_f.write("notes\tSummary generated from associationByOverallIndirect JSON files\n")

print("Done.")
print("Rows processed:", total_rows)
print("Unique genes:", len(gene_stats))
print("Bad lines:", bad_lines)
print("Gene summary:", SUMMARY_TSV)
print("General summary:", GENERAL_TXT)
print("Metadata:", METADATA_TXT)
