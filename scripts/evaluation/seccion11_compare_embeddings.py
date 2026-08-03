#!/usr/bin/env python3
"""
Seccion 11 (TFM graph_v4): comparacion GEOMETRICA (no coordenada a
coordenada) entre los 9 embeddings aprendidos del catalogo unificado
(Secciones 9-10) y los embeddings externos de proteinas UniProt-ProtT5
/ GenePT-Model3 (Seccion 8).

Por que geometrica y no directa: los espacios tienen dimensiones
distintas (128 vs 1.024/3.072) y pueden estar rotados/reflejados/
escalados sin perder equivalencia relacional -- comparar coordenada a
coordenada (p.ej. correlacion por dimension) no tiene sentido. En su
lugar se compara la ESTRUCTURA RELACIONAL entre genes: que tan parecidas
son las relaciones de similitud/vecindad entre genes en un espacio
frente al otro.

CONJUNTO COMUN DE GENES (regla explicita del protocolo): para cada
comparacion (modelo, externo), el conjunto base es la interseccion de
node_id's presentes en ambos. CKA y Ridge (que no requieren muestreo de
pares, son baratos sobre matrices n x d) usan ese conjunto COMPLETO.
Spearman (correlacion de matrices de similitud por pares) y el
solapamiento de vecinos (Jaccard) SI requieren acotar el numero de
pares/consultas por costo computacional -- para ambas se extrae UNA
UNICA submuestra reproducible (seed=42) del conjunto comun, y esa MISMA
submuestra se usa para las dos metricas (tal como exige el protocolo:
"usar exactamente el mismo conjunto de genes para todas las metricas de
una comparacion"). Esta eleccion se documenta explicitamente aqui y en
el entregable final, no es un desvio silencioso.

Metricas (nunca combinadas en un indice unico):
  1. Spearman de matrices de similitud coseno por pares (submuestra).
  2. Jaccard de vecinos mas cercanos, k=10,25,50 (mismas consultas que 1,
     pero buscando vecinos sobre el conjunto comun COMPLETO).
  3. CKA lineal (Kornblith et al. 2019), formulacion eficiente sin
     materializar matrices n x n.
  4. Regresion Ridge multisalida con validacion cruzada de 5 folds sobre
     genes, prediciendo el embedding EXTERNO a partir del APRENDIDO;
     R^2 promedio en folds de test. Presentado como correspondencia
     lineal, no como equivalencia biologica (ver entregable final).

Uso (en un nodo de computo, nunca en el login node):
    python3 seccion11_compare_embeddings.py \
        --out-dir data/processed/analysis/seccion11_comparacion_embeddings
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

SAMPLE_SIZE_PAIRWISE = 1500  # genes, para Spearman de pares y consultas kNN
K_VALUES = [10, 25, 50]
SEED = 42
MIN_COMMON_GENES = 30  # por debajo de esto, comparacion omitida (no se fuerza)

CATALOG_DIR = "data/processed/analysis/embedding_catalog_v4"
EXTERNAL_DIR = "data/processed/analysis/external_embeddings_v4"

MODELS = [
    ("spectral", f"{CATALOG_DIR}/spectral_structural_seed42.tsv"),
    ("node2vec", f"{CATALOG_DIR}/node2vec_structural_seed42.tsv"),
    ("hin2vec", f"{CATALOG_DIR}/hin2vec_structural_seed42.tsv"),
    ("gae", f"{CATALOG_DIR}/gae_structural_seed42.tsv"),
    ("gat", f"{CATALOG_DIR}/gat_structural_seed42.tsv"),
    ("gae_uniprot_fixed", f"{CATALOG_DIR}/gae_uniprot_fixed_seed42.tsv"),
    ("gae_uniprot_finetuned", f"{CATALOG_DIR}/gae_uniprot_finetuned_seed42.tsv"),
    ("gat_uniprot_fixed", f"{CATALOG_DIR}/gat_uniprot_fixed_seed42.tsv"),
    ("gat_uniprot_finetuned", f"{CATALOG_DIR}/gat_uniprot_finetuned_seed42.tsv"),
]

EXTERNALS = [
    ("UniProt_ProtT5", f"{EXTERNAL_DIR}/uniprot_embeddings_subgraph.tsv"),
    ("GenePT_Model3", f"{EXTERNAL_DIR}/genept_embeddings_subgraph.tsv"),
]


def log(msg):
    print(f"[seccion11] {msg}", flush=True)


def load_embeddings(path, is_external):
    df = pd.read_csv(path, sep="\t")
    e_cols = sorted((c for c in df.columns if c.startswith("e_")), key=lambda c: int(c.split("_")[1]))
    if not is_external:
        df = df[df["node_type"] == "GEN"]
    ids = df["node_id"].to_numpy()
    mat = df[e_cols].to_numpy(dtype=np.float32)
    return dict(zip(ids, mat)), len(e_cols)


def l2_normalize(mat):
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return mat / norms


def spearman_pairwise_similarity(A_norm_sample, B_norm_sample):
    sim_A = A_norm_sample @ A_norm_sample.T
    sim_B = B_norm_sample @ B_norm_sample.T
    iu = np.triu_indices(sim_A.shape[0], k=1)
    rho, pval = spearmanr(sim_A[iu], sim_B[iu])
    return float(rho), float(pval), int(len(iu[0]))


def jaccard_neighbor_overlap(sim_own, sim_ext, k):
    """sim_own, sim_ext: (m_consultas, n_candidatos), auto-similitud ya
    anulada (-inf) donde corresponda. Devuelve la media del indice de
    Jaccard de los k vecinos mas cercanos por consulta."""
    m = sim_own.shape[0]
    jaccards = np.empty(m)
    for i in range(m):
        own_top = set(np.argpartition(-sim_own[i], k)[:k].tolist())
        ext_top = set(np.argpartition(-sim_ext[i], k)[:k].tolist())
        inter = len(own_top & ext_top)
        union = len(own_top | ext_top)
        jaccards[i] = inter / union if union else 0.0
    return float(jaccards.mean())


def linear_cka(X, Y):
    """CKA lineal (Kornblith, Norouzi, Lee, Hinton 2019), formulacion
    eficiente que evita materializar matrices n x n:
    CKA = ||Yc^T Xc||_F^2 / (||Xc^T Xc||_F * ||Yc^T Yc||_F)
    con Xc, Yc centradas por columnas."""
    Xc = X - X.mean(axis=0, keepdims=True)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    xty = Xc.T @ Yc
    xtx = Xc.T @ Xc
    yty = Yc.T @ Yc
    num = np.linalg.norm(xty, ord="fro") ** 2
    den = np.linalg.norm(xtx, ord="fro") * np.linalg.norm(yty, ord="fro")
    return float(num / den) if den > 0 else float("nan")


def ridge_cv_r2(X, Y, seed=SEED, n_splits=5, alpha=1.0):
    """Ridge multisalida X(aprendido) -> Y(externo), 5-fold CV sobre
    genes. R^2 (uniform_average sobre dimensiones de salida, default de
    sklearn) promediado entre folds."""
    n_splits_eff = min(n_splits, len(X))
    if n_splits_eff < 2:
        return float("nan"), float("nan")
    kf = KFold(n_splits=n_splits_eff, shuffle=True, random_state=seed)
    scores = []
    for train_idx, test_idx in kf.split(X):
        model = Ridge(alpha=alpha)
        model.fit(X[train_idx], Y[train_idx])
        scores.append(model.score(X[test_idx], Y[test_idx]))
    return float(np.mean(scores)), float(np.std(scores))


def compare_one(model_name, model_emb, ext_name, ext_emb, out_rows, pca_examples):
    common = sorted(set(model_emb.keys()) & set(ext_emb.keys()))
    n_common = len(common)
    log(f"{model_name} vs {ext_name}: {n_common:,} genes comunes")
    if n_common < MIN_COMMON_GENES:
        out_rows.append(dict(
            modelo=model_name, externo=ext_name, n_genes_comunes=n_common,
            n_genes_muestra_pares=None, n_pares_spearman=None,
            spearman=None, spearman_pval=None,
            jaccard_k10=None, jaccard_k25=None, jaccard_k50=None,
            cka_lineal=None, ridge_r2=None, ridge_r2_std=None,
            nota=f"menos de {MIN_COMMON_GENES} genes comunes, comparacion omitida",
        ))
        return

    M = np.stack([model_emb[g] for g in common])
    E = np.stack([ext_emb[g] for g in common])

    # --- CKA y Ridge: conjunto COMPLETO de genes comunes ---
    cka = linear_cka(M, E)
    ridge_r2_mean, ridge_r2_std = ridge_cv_r2(M, E, seed=SEED)

    # --- Spearman y Jaccard: UNA submuestra reproducible, misma para ambas ---
    rng = np.random.default_rng(SEED)
    sample_size = min(SAMPLE_SIZE_PAIRWISE, n_common)
    sample_idx = rng.choice(n_common, size=sample_size, replace=False)

    M_norm = l2_normalize(M)
    E_norm = l2_normalize(E)
    M_sample = M_norm[sample_idx]
    E_sample = E_norm[sample_idx]

    rho, pval, n_pairs = spearman_pairwise_similarity(M_sample, E_sample)

    sim_own = M_sample @ M_norm.T   # (sample, n_common)
    sim_ext = E_sample @ E_norm.T
    for j, gi in enumerate(sample_idx):
        sim_own[j, gi] = -np.inf
        sim_ext[j, gi] = -np.inf

    jacc = {k: jaccard_neighbor_overlap(sim_own, sim_ext, k) for k in K_VALUES}

    out_rows.append(dict(
        modelo=model_name, externo=ext_name, n_genes_comunes=n_common,
        n_genes_muestra_pares=sample_size, n_pares_spearman=n_pairs,
        spearman=rho, spearman_pval=pval,
        jaccard_k10=jacc[10], jaccard_k25=jacc[25], jaccard_k50=jacc[50],
        cka_lineal=cka, ridge_r2=ridge_r2_mean, ridge_r2_std=ridge_r2_std, nota="",
    ))
    pca_examples[(model_name, ext_name)] = (M, E, common)


def make_heatmap(df, value_col, title, out_path, fmt="{:.3f}"):
    pivot = df.pivot(index="modelo", columns="externo", values=value_col)
    pivot = pivot.reindex([m for m, _ in MODELS])
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    im = ax.imshow(pivot.values.astype(float), cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=20, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if pd.notna(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center",
                        color="white" if v < np.nanmax(pivot.values) * 0.6 else "black", fontsize=8)
    ax.set_title(title, fontsize=11, wrap=True)
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def make_jaccard_bars(df, out_path):
    fig, axes = plt.subplots(1, len(EXTERNALS), figsize=(13, 5), sharey=True)
    x = np.arange(len(MODELS))
    width = 0.25
    for ax, (ext_name, _) in zip(axes, EXTERNALS):
        sub = df[df["externo"] == ext_name].set_index("modelo").reindex([m for m, _ in MODELS])
        for i, k in enumerate(K_VALUES):
            ax.bar(x + (i - 1) * width, sub[f"jaccard_k{k}"].astype(float), width, label=f"k={k}")
        ax.set_xticks(x)
        ax.set_xticklabels([m for m, _ in MODELS], rotation=45, ha="right")
        ax.set_title(f"Solapamiento de vecinos (Jaccard) vs {ext_name}")
        ax.set_ylabel("Índice de Jaccard (media)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def make_cka_plot(df, out_path):
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(MODELS))
    width = 0.35
    for i, (ext_name, _) in enumerate(EXTERNALS):
        sub = df[df["externo"] == ext_name].set_index("modelo").reindex([m for m, _ in MODELS])
        ax.bar(x + (i - 0.5) * width, sub["cka_lineal"].astype(float), width, label=ext_name)
    ax.set_xticks(x)
    ax.set_xticklabels([m for m, _ in MODELS], rotation=45, ha="right")
    ax.set_ylabel("CKA lineal")
    ax.set_title("CKA lineal: modelo aprendido vs UniProt / GenePT")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def make_pca_exploratory(pca_examples, choice_key, out_path, title):
    if choice_key not in pca_examples:
        log(f"AVISO: no hay datos para PCA exploratorio de {choice_key}, se omite figura.")
        return
    M, E, common = pca_examples[choice_key]
    n_show = min(2000, len(common))
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(common), size=n_show, replace=False)
    pca_m = PCA(n_components=2, random_state=SEED).fit_transform(M[idx])
    pca_e = PCA(n_components=2, random_state=SEED).fit_transform(E[idx])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].scatter(pca_m[:, 0], pca_m[:, 1], s=4, alpha=0.5)
    axes[0].set_title(f"PCA — {choice_key[0]} (aprendido)")
    axes[1].scatter(pca_e[:, 0], pca_e[:, 1], s=4, alpha=0.5, color="darkorange")
    axes[1].set_title(f"PCA — {choice_key[1]} (externo)")
    fig.suptitle(title + "\n(exploratorio -- no es evidencia de equivalencia entre espacios)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    log("Cargando embeddings externos...")
    externals = {}
    for name, path in EXTERNALS:
        emb, dim = load_embeddings(path, is_external=True)
        externals[name] = emb
        log(f"  {name}: {len(emb):,} genes, dim={dim}")

    out_rows = []
    pca_examples = {}
    for model_name, path in MODELS:
        log(f"Cargando {model_name}...")
        model_emb, dim = load_embeddings(path, is_external=False)
        log(f"  {model_name}: {len(model_emb):,} nodos GEN, dim={dim}")
        for ext_name, ext_emb in externals.items():
            compare_one(model_name, model_emb, ext_name, ext_emb, out_rows, pca_examples)

    out_df = pd.DataFrame(out_rows)
    metrics_path = f"{args.out_dir}/seccion11_metricas_comparacion.tsv"
    out_df.to_csv(metrics_path, sep="\t", index=False)
    log(f"Tabla de metricas -> {metrics_path}")

    # --- Figuras ---
    valid_df = out_df.dropna(subset=["spearman"])
    if not valid_df.empty:
        make_heatmap(valid_df, "spearman", "Spearman (similitud por pares): aprendido vs externo",
                     f"{args.out_dir}/seccion11_heatmap_spearman.png")
        make_jaccard_bars(valid_df, f"{args.out_dir}/seccion11_barras_jaccard.png")
        make_cka_plot(valid_df, f"{args.out_dir}/seccion11_cka.png")

        # PCA exploratorio: mejor y peor CKA vs UniProt (solo ilustrativo)
        uniprot_rows = valid_df[valid_df["externo"] == "UniProt_ProtT5"].dropna(subset=["cka_lineal"])
        if not uniprot_rows.empty:
            best = uniprot_rows.loc[uniprot_rows["cka_lineal"].idxmax(), "modelo"]
            worst = uniprot_rows.loc[uniprot_rows["cka_lineal"].idxmin(), "modelo"]
            make_pca_exploratory(pca_examples, (best, "UniProt_ProtT5"),
                                 f"{args.out_dir}/seccion11_pca_mejor_cka_uniprot.png",
                                 f"Mejor CKA vs UniProt: {best}")
            make_pca_exploratory(pca_examples, (worst, "UniProt_ProtT5"),
                                 f"{args.out_dir}/seccion11_pca_peor_cka_uniprot.png",
                                 f"Peor CKA vs UniProt: {worst}")
    else:
        log("AVISO: ninguna comparacion tuvo suficientes genes comunes -- no se generan figuras.")

    log("Listo.")


if __name__ == "__main__":
    main()
