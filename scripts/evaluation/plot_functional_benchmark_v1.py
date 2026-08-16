#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')  # Headless mode fix to prevent display crashes

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path("data/processed/analysis/functional_benchmark_v1")
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

ORDER = [
    "genept_raw", "node2vec", "uniprot_raw",
    "gat_uniprot_fixed", "gae_uniprot_finetuned",
    "gat_uniprot_finetuned", "gae_uniprot_fixed",
    "spectral", "gae_structural", "gat_structural", "hin2vec",
]
LABELS = {
    "genept_raw": "GenePT",
    "node2vec": "Node2Vec",
    "uniprot_raw": "UniProt",
    "gat_uniprot_fixed": "GAT + UniProt\nfixed",
    "gae_uniprot_finetuned": "GAE + UniProt\nfine-tuned",
    "gat_uniprot_finetuned": "GAT + UniProt\nfine-tuned",
    "gae_uniprot_fixed": "GAE + UniProt\nfixed",
    "spectral": "Spectral",
    "gae_structural": "GAE\nstructural",
    "gat_structural": "GAT\nstructural",
    "hin2vec": "HIN2Vec",
}
PALETTE = {
    "genept_raw": "#7b3294", "node2vec": "#008837", "uniprot_raw": "#2166ac",
    "gat_uniprot_fixed": "#c2a5cf", "gae_uniprot_finetuned": "#f4a582",
    "gat_uniprot_finetuned": "#d6604d", "gae_uniprot_fixed": "#fddbc7",
    "spectral": "#636363", "gae_structural": "#9ecae1",
    "gat_structural": "#6baed6", "hin2vec": "#bdbdbd",
}

sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)

go = pd.read_csv(ROOT / "go_all_term_results.tsv", sep="\t")
go = go[go["status"] == "ok"].copy()
go["label"] = go["source"].map(LABELS)
go["source"] = pd.Categorical(go["source"], categories=ORDER, ordered=True)
go = go.sort_values("source")

# Figure 1: GO distribution
fig, ax = plt.subplots(figsize=(12, 6.6))
sns.boxplot(
    data=go, x="source", y="test_auprc", order=ORDER,
    hue="source", palette=PALETTE, legend=False, dodge=False,
    width=0.68, showfliers=False, linewidth=1, ax=ax,
)
sns.stripplot(
    data=go, x="source", y="test_auprc", order=ORDER,
    color="#222222", alpha=0.28, size=2.6, jitter=0.22, ax=ax,
)
means = go.groupby("source", observed=True)["test_auprc"].mean()
for i, source in enumerate(ORDER):
    ax.scatter(i, means[source], color="white", edgecolor="black",
               linewidth=0.8, s=38, zorder=4)
ax.set_xlabel("")
ax.set_ylabel("AUPRC en holdout")
ax.set_title("Predicción funcional GO-BP: distribución entre 56 términos")
ax.set_xticks(range(len(ORDER)), [LABELS[s] for s in ORDER],
              rotation=42, ha="right")
ax.set_ylim(0, 1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_go_auprc_distribution.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "figure_go_auprc_distribution.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: BioGRID grouped bars
bg = pd.read_csv(ROOT / "biogrid_all_fold_results.tsv", sep="\t")
bg = bg[bg["status"] == "ok"].copy()
summary = (
    bg.groupby(["source", "task"], as_index=False)["test_auprc"]
    .mean()
)
summary["source"] = pd.Categorical(summary["source"], categories=ORDER, ordered=True)
summary = summary.sort_values("source")
summary["task"] = summary["task"].map({
    "synthetic_lethality": "Letalidad sintética",
    "negative_genetic": "Interacción genética negativa",
})

fig, ax = plt.subplots(figsize=(12, 6.6))
sns.barplot(
    data=summary, x="source", y="test_auprc", hue="task",
    order=ORDER,
    hue_order=["Letalidad sintética", "Interacción genética negativa"],
    palette=["#b2182b", "#2166ac"], errorbar=None, ax=ax,
)
ax.set_xlabel("")
ax.set_ylabel("AUPRC media en test")
ax.set_title("Predicción de interacciones genéticas externas de BioGRID")
ax.set_xticks(range(len(ORDER)), [LABELS[s] for s in ORDER],
              rotation=42, ha="right")
ax.legend(title="")
ax.set_ylim(0, 0.55)
fig.tight_layout()
fig.savefig(OUT / "figure_biogrid_auprc.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "figure_biogrid_auprc.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 3: unified heatmap
master = pd.read_csv(ROOT / "functional_benchmark_master_summary.tsv", sep="\t")
master = master.set_index("source").loc[ORDER]
heat = master[[
    "go_auprc_mean",
    "synthetic_lethality_test_auprc",
    "negative_genetic_test_auprc",
]].copy()
heat.columns = ["GO-BP", "BioGRID SL", "BioGRID NG"]
heat.index = [LABELS[s] for s in heat.index]

fig, ax = plt.subplots(figsize=(7.4, 7.2))
sns.heatmap(
    heat, annot=True, fmt=".3f", cmap="YlGnBu", vmin=0, vmax=0.8,
    linewidths=0.6, linecolor="white",
    cbar_kws={"label": "AUPRC"}, ax=ax,
)
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("Comparación de rendimiento funcional downstream")
fig.tight_layout()
fig.savefig(OUT / "figure_functional_auprc_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "figure_functional_auprc_heatmap.pdf", bbox_inches="tight")
plt.close(fig)

print(f"Saved figures to: {OUT}")
