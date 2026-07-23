#!/usr/bin/env python
"""
Identifica y flaggea ambiguedades en el mapeo ENSG <-> UniProt:
  - ENSG que mapean a multiples UniProt accessions (isoformas)
  - UniProt accessions compartidos por multiples ENSG (paralogos, genes
    en region pseudoautosomica, etc.)

Toma como input el archivo ya generado por check_uniprot_coverage.py
(ensg_uniprot_mapping_full.tsv) - no recalcula el cruce desde cero.

Uso:
    python flag_ambiguous_mappings.py \
        --mapping-full data/processed/analysis/uniprot_coverage/ensg_uniprot_mapping_full.tsv \
        --out-dir data/processed/analysis/uniprot_coverage

Ejecutar en nodo de computo (srun/sbatch), no en login node, por uso de pandas.
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mapping-full", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.mapping_full, sep="\t")

    mapped = df[df["uniprot"].notna()].copy()

    n_uniprot_per_ensg = mapped.groupby("ensg")["uniprot"].nunique()
    n_ensg_per_uniprot = mapped.groupby("uniprot")["ensg"].nunique()

    mapped["n_uniprot_para_este_ensg"] = mapped["ensg"].map(n_uniprot_per_ensg)
    mapped["n_ensg_para_este_uniprot"] = mapped["uniprot"].map(n_ensg_per_uniprot)

    mapped["flag_ensg_multi_uniprot"] = mapped["n_uniprot_para_este_ensg"] > 1
    mapped["flag_uniprot_multi_ensg"] = mapped["n_ensg_para_este_uniprot"] > 1

    # Reincorporar filas sin mapeo para no perder genes del total
    unmapped = df[df["uniprot"].isna()].copy()
    for col in [
        "n_uniprot_para_este_ensg",
        "n_ensg_para_este_uniprot",
        "flag_ensg_multi_uniprot",
        "flag_uniprot_multi_ensg",
    ]:
        unmapped[col] = pd.NA

    flagged_df = pd.concat([mapped, unmapped], ignore_index=True)

    flagged_path = out_dir / "ensg_uniprot_mapping_flagged.tsv"
    flagged_df.to_csv(flagged_path, sep="\t", index=False)

    # Listados separados para revision humana
    ensg_multi = (
        mapped[mapped["flag_ensg_multi_uniprot"]]
        .groupby("ensg")["uniprot"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .reset_index()
        .rename(columns={"uniprot": "uniprot_accessions"})
    )
    ensg_multi["n_uniprot"] = ensg_multi["uniprot_accessions"].str.count(";") + 1

    uniprot_multi = (
        mapped[mapped["flag_uniprot_multi_ensg"]]
        .groupby("uniprot")["ensg"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .reset_index()
        .rename(columns={"ensg": "ensg_genes"})
    )
    uniprot_multi["n_ensg"] = uniprot_multi["ensg_genes"].str.count(";") + 1

    ensg_multi_path = out_dir / "ensg_con_multiples_uniprot.tsv"
    uniprot_multi_path = out_dir / "uniprot_compartido_por_multiples_ensg.tsv"

    ensg_multi.to_csv(ensg_multi_path, sep="\t", index=False)
    uniprot_multi.to_csv(uniprot_multi_path, sep="\t", index=False)

    n_ensg_afectados_por_compartir = mapped.loc[
        mapped["flag_uniprot_multi_ensg"], "ensg"
    ].nunique()

    print(f"ENSG con multiples UniProt (isoformas): {len(ensg_multi)}")
    print(f"UniProt compartidos por multiples ENSG: {len(uniprot_multi)}")
    print(f"Genes ENSG afectados por compartir UniProt con otro: {n_ensg_afectados_por_compartir}")
    print()
    print(f"Archivo flaggeado completo: {flagged_path}")
    print(f"Listado ENSG con multi-UniProt: {ensg_multi_path}")
    print(f"Listado UniProt compartido: {uniprot_multi_path}")


if __name__ == "__main__":
    main()
