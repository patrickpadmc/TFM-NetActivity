#!/usr/bin/env python3
"""
Codigo de evaluacion comun para todas las tareas y metodos de embedding
del TFM (graph_v4). Se usa igual para Spectral Embedding, Node2Vec,
HIN2Vec, GAE y GAT-AE, de forma que la comparacion entre metodos no
dependa de diferencias de implementacion en la parte de evaluacion.

Piezas:
  - build_features(): construye el vector de features de un par de nodos
    a partir de sus embeddings, mediante producto de Hadamard + diferencia
    absoluta (concatenados), el minimo exigido por el protocolo.
  - train_classifier(): regresion logistica con ajuste de hiperparametro
    (C) por AUPRC en validacion.
  - evaluate(): AUPRC (metrica primaria) y AUROC (metrica secundaria).
  - ranked_evaluation(): Hits@k / MRR mediante ranking contra un pool de
    negativos muestreados (protocolo de "ranking contra N negativos",
    necesario porque rankear cada positivo de test contra el universo
    completo de candidatos no es computacionalmente viable a esta escala).

IMPORTANTE (metrica primaria): AUPRC/Average Precision. AUROC es
secundaria. La perdida de reconstruccion de GAE/GAT-AE NUNCA debe
presentarse como metrica de comparacion principal -- solo como
diagnostico interno de esos dos metodos.

Decisiones de diseno que requieren validacion del usuario antes de
usarse en resultados finales (ver docstring de ranked_evaluation):
  - tamano del pool de negativos de ranking (por defecto 99, protocolo
    habitual en la literatura de "1 positivo vs N negativos").
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

DEFAULT_C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)


def build_features(embeddings: dict, pairs, missing_policy: str = "zero"):
    """
    embeddings: dict nodo_id -> np.ndarray (misma dimension para todos los
        nodos de un mismo metodo/ejecucion de embedding).
    pairs: iterable de tuplas (n1, n2).
    missing_policy: que hacer si un nodo no tiene embedding
        ("zero" -> vector de ceros, "skip" -> se omite el par).

    Devuelve (X, kept_mask) donde X es un array (n_pairs_validos, 2*dim)
    con el producto de Hadamard concatenado con la diferencia absoluta,
    y kept_mask es un array booleano indicando que pares de la entrada
    original se conservaron (relevante solo si missing_policy="skip").
    """
    pairs = list(pairs)
    if not pairs:
        return np.empty((0, 0)), np.array([], dtype=bool)

    dim = len(next(iter(embeddings.values())))
    feats = []
    kept = []
    for n1, n2 in pairs:
        e1 = embeddings.get(n1)
        e2 = embeddings.get(n2)
        if e1 is None or e2 is None:
            if missing_policy == "skip":
                kept.append(False)
                continue
            e1 = np.zeros(dim) if e1 is None else e1
            e2 = np.zeros(dim) if e2 is None else e2
        hadamard = e1 * e2
        abs_diff = np.abs(e1 - e2)
        feats.append(np.concatenate([hadamard, abs_diff]))
        kept.append(True)
    X = np.vstack(feats) if feats else np.empty((0, 2 * dim))
    return X, np.array(kept, dtype=bool)


def train_classifier(X_train, y_train, X_val, y_val, c_grid=DEFAULT_C_GRID, random_state=42):
    """
    Ajusta una regresion logistica probando cada valor de c_grid,
    seleccionando el que maximiza AUPRC en validacion. Devuelve el
    modelo reentrenado sobre train+val con el mejor C, junto con un
    dict de diagnostico (mejor C, AUPRC de validacion por candidato).
    """
    best_c, best_auprc, val_scores = None, -1.0, {}
    for c in c_grid:
        clf = LogisticRegression(C=c, max_iter=2000, random_state=random_state)
        clf.fit(X_train, y_train)
        val_pred = clf.predict_proba(X_val)[:, 1]
        auprc = average_precision_score(y_val, val_pred)
        val_scores[c] = auprc
        if auprc > best_auprc:
            best_auprc, best_c = auprc, c

    X_full = np.vstack([X_train, X_val])
    y_full = np.concatenate([y_train, y_val])
    final_clf = LogisticRegression(C=best_c, max_iter=2000, random_state=random_state)
    final_clf.fit(X_full, y_full)

    diagnostics = {"best_c": best_c, "val_auprc_por_c": val_scores, "val_auprc_mejor_c": best_auprc}
    return final_clf, diagnostics


def evaluate(model, X_test, y_test):
    """AUPRC (primaria) y AUROC (secundaria) sobre el conjunto de test."""
    scores = model.predict_proba(X_test)[:, 1]
    return {
        "auprc": average_precision_score(y_test, scores),
        "auroc": roc_auc_score(y_test, scores),
        "n_test": int(len(y_test)),
        "n_pos_test": int(y_test.sum()),
    }


def ranked_evaluation(
    model,
    embeddings: dict,
    test_positive_pairs,
    node_pool_t2,
    rng,
    k: int = 10,
    n_negatives_per_query: int = 99,
    forbidden_pairs: set | None = None,
):
    """
    Hits@k y MRR mediante ranking de cada positivo de test contra un
    pool de n_negatives_per_query candidatos muestreados aleatoriamente
    del universo de nodos de tipo t2 (protocolo estandar "1 positivo
    vs N negativos", NO ranking contra el universo completo -- a esta
    escala (hasta 26,330 x 23,697 pares posibles) el ranking exhaustivo
    por consulta es computacionalmente inviable).

    forbidden_pairs: conjunto de pares (n1, n2) conocidos como positivos
    reales (train+val+test) para ese universo t1-t2, de forma que un
    negativo muestreado nunca sea en realidad un enlace verdadero.

    Devuelve dict con hits_at_k, mrr, n_queries_evaluadas.

    NOTA DE DISENO A VALIDAR: n_negatives_per_query=99 es el valor por
    defecto; debe confirmarse con el usuario antes de reportarse como
    resultado final del TFM, ya que el tamano del pool afecta la escala
    absoluta de Hits@k (no a AUPRC/AUROC, que son las metricas primaria
    y secundaria).
    """
    forbidden_pairs = forbidden_pairs or set()
    ranks = []
    for n1, n2_true in test_positive_pairs:
        if n1 not in embeddings or n2_true not in embeddings:
            continue
        candidates = [n2_true]
        guard = 0
        while len(candidates) < n_negatives_per_query + 1 and guard < 20:
            guard += 1
            sample = rng.choice(node_pool_t2, size=n_negatives_per_query * 2)
            for cand in sample.tolist():
                if len(candidates) >= n_negatives_per_query + 1:
                    break
                if cand == n2_true or (n1, cand) in forbidden_pairs:
                    continue
                if cand in embeddings:
                    candidates.append(cand)
        if len(candidates) < 2:
            continue

        e1 = embeddings[n1]
        feats = np.vstack([
            np.concatenate([e1 * embeddings[c], np.abs(e1 - embeddings[c])])
            for c in candidates
        ])
        scores = model.predict_proba(feats)[:, 1]
        order = np.argsort(-scores)
        rank = int(np.where(order == 0)[0][0]) + 1
        ranks.append(rank)

    if not ranks:
        return {"hits_at_k": None, "mrr": None, "n_queries": 0, "k": k}

    ranks = np.array(ranks)
    hits_at_k = float((ranks <= k).mean())
    mrr = float((1.0 / ranks).mean())
    return {"hits_at_k": hits_at_k, "mrr": mrr, "n_queries": int(len(ranks)), "k": k}
