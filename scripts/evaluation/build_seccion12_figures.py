#!/usr/bin/env python3
"""Figuras de Seccion 12 (TFM graph_v4):
  1. Matriz visual de AUPRC (fuente x relacion).
  2. Curvas Precision-Recall de las condiciones principales sobre
     GEN-ppi-GEN (la unica relacion comun a las 11 fuentes, incluidas
     UniProt/GenePT crudos).
  3. Curvas ROC, misma seleccion.
  4. Curvas PR de las 9 fuentes de grafo completo sobre GEN-ass-DIS (la
     tarea de mayor volumen), para comparar estructurales vs UniProt-init.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = "data/processed/analysis/seccion12_svm_link_prediction"

MODELS_ALL = ["spectral", "node2vec", "hin2vec", "gae", "gat",
              "gae_uniprot_fixed", "gae_uniprot_finetuned",
              "gat_uniprot_fixed", "gat_uniprot_finetuned"]
MODELS_WITH_RAW = MODELS_ALL + ["uniprot_raw", "genept_raw"]


def load_curve(source, relation):
    return np.load(f"{OUT_DIR}/{source}__{relation}__curves.npz")


def plot_auprc_matrix():
    pivot = pd.read_csv(f"{OUT_DIR}/seccion12_matriz_auprc.tsv", sep="\t", index_col=0)
    pivot = pivot.reindex(MODELS_WITH_RAW)
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(pivot.values.astype(float), cmap="viridis", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v < 0.6 else "black", fontsize=8)
    ax.set_title("AUPRC (test) por fuente de embedding y relación — SVM (Sección 12)", fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.75, label="AUPRC")
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/seccion12_matriz_auprc.png", dpi=150)
    plt.close(fig)


def plot_pr_roc(relation, sources, out_prefix, title_suffix):
    fig_pr, ax_pr = plt.subplots(figsize=(7, 6))
    fig_roc, ax_roc = plt.subplots(figsize=(7, 6))
    cmap = plt.get_cmap("tab20")
    for i, src in enumerate(sources):
        try:
            c = load_curve(src, relation)
        except FileNotFoundError:
            continue
        color = cmap(i % 20)
        ax_pr.plot(c["recall"], c["precision"], label=src, color=color, linewidth=1.3)
        ax_roc.plot(c["fpr"], c["tpr"], label=src, color=color, linewidth=1.3)
    ax_pr.set_xlabel("Recall"); ax_pr.set_ylabel("Precision")
    ax_pr.set_title(f"Precision-Recall — {relation} {title_suffix}", fontsize=10)
    ax_pr.legend(fontsize=7, loc="lower left")
    fig_pr.tight_layout()
    fig_pr.savefig(f"{OUT_DIR}/{out_prefix}_pr.png", dpi=150)
    plt.close(fig_pr)

    ax_roc.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5)
    ax_roc.set_xlabel("FPR (1 - especificidad)"); ax_roc.set_ylabel("TPR (recall)")
    ax_roc.set_title(f"ROC — {relation} {title_suffix}", fontsize=10)
    ax_roc.legend(fontsize=7, loc="lower right")
    fig_roc.tight_layout()
    fig_roc.savefig(f"{OUT_DIR}/{out_prefix}_roc.png", dpi=150)
    plt.close(fig_roc)


def main():
    plot_auprc_matrix()
    print("Matriz AUPRC -> seccion12_matriz_auprc.png")

    plot_pr_roc("GEN-ppi-GEN", MODELS_WITH_RAW, "seccion12_gen_ppi_gen",
                "(las 11 fuentes, unica relacion comun a UniProt/GenePT crudos)")
    print("PR/ROC GEN-ppi-GEN -> seccion12_gen_ppi_gen_{pr,roc}.png")

    plot_pr_roc("GEN-ass-DIS", MODELS_ALL, "seccion12_gen_ass_dis",
                "(9 fuentes de grafo completo, tarea de mayor volumen)")
    print("PR/ROC GEN-ass-DIS -> seccion12_gen_ass_dis_{pr,roc}.png")

    print("Listo.")


if __name__ == "__main__":
    main()
