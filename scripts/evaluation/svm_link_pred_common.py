#!/usr/bin/env python3
"""
Seccion 12 (TFM graph_v4): funciones comunes para el pipeline de
prediccion de enlaces con SVM, por (fuente de embedding, relacion).

Construccion de vectores de par (protocolo Seccion 12: "Producto de
Hadamard. Diferencia absoluta. Concatenacion, si es viable"):
  - Hadamard (e1 * e2) y diferencia absoluta (|e1 - e2|) SIEMPRE se
    incluyen (dimension = dim del embedding en cada caso).
  - Concatenacion cruda [e1, e2] se incluye SOLO cuando dim del
    embedding es baja (dim <= CONCAT_DIM_THRESHOLD): para los
    embeddings crudos UniProt (1024-d) / GenePT (3072-d), anadir la
    concatenacion multiplicaria el tamano del vector de features por 4
    en vez de por 2 (p.ej. GenePT: 4*3072=12.288 features x ~400.000
    pares de GEN-ppi-GEN séria inviable en memoria). Esta es la lectura
    operativa explicita de "si es viable" del protocolo.
  - missing_policy="skip": si falta el embedding de alguno de los dos
    nodos del par, el par se descarta (nunca se inventa un embedding
    externo). Se cuenta y reporta cuantos pares se descartaron por
    fuente/relacion -- relevante sobre todo para UniProt/GenePT crudos,
    que no cubren el 100% de los nodos GEN.

Escala del SVM: los conjuntos de entrenamiento de esta seccion llegan
hasta ~440.000 pares (GEN-ass-TIS, GEN-ppi-GEN). Un SVM de kernel RBF
EXACTO (sklearn.svm.SVC) tiene coste O(n^2)-O(n^3) y es inviable a esta
escala. En su lugar se usa una aproximacion de kernel RBF mediante
mapeo de Nystroem + SVM lineal (SGDClassifier, ver nota de escala mas abajo) sobre las features
transformadas -- coste lineal en el numero de muestras, practica
estandar para "SVM de kernel a gran escala". Esto se documenta aqui
explicitamente: NO es un SVM RBF exacto.

Seleccion lineal vs RBF-aproximado: se entrenan ambos (grid de C para
el lineal; grid de gamma x C para el RBF-aproximado), se comparan por
AUPRC de VALIDACION, y se conserva el que tenga mejor AUPRC de
validacion -- asi se cumple "si aporta mejora verificable en
validacion, probar RBF" sin decidir de antemano si merece la pena
probarlo.

Umbral: se elige en VALIDACION maximizando F1, nunca en test. AUPRC y
AUROC no requieren umbral (se calculan sobre el score continuo,
model.decision_function).

NOTA DE ESCALA (anadida tras el primer intento real en cluster): con
sklearn.svm.LinearSVC (liblinear) el ajuste NO convergio en un tiempo
razonable sobre las tareas mas grandes (~400.000 pares de train,
partition=short, 4h) -- ni siquiera termino el primer valor de C de la
primera tarea. liblinear escala mal en este regimen. Se sustituye por
sklearn.linear_model.SGDClassifier(loss="hinge") -- SVM lineal
optimizado por descenso de gradiente estocastico con parada temprana,
la alternativa estandar de sklearn para SVM lineal a gran escala de
muestras (ver documentacion de sklearn sobre eleccion de solver segun
tamano de dataset). alpha se calcula como 1/(C * n_train), la
correspondencia aproximada estandar entre el hiperparametro C de un SVM
y alpha de SGDClassifier con perdida hinge.
"""
import time

import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.kernel_approximation import Nystroem
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    f1_score,
    precision_score,
    recall_score,
    precision_recall_curve,
    roc_curve,
)

CONCAT_DIM_THRESHOLD = 256  # concatenacion cruda [e1,e2] solo si dim <= esto
C_GRID = (0.01, 0.1, 1.0, 10.0)
GAMMA_GRID = (None, 0.1, 1.0)  # None = default de Nystroem (~1/n_features)
N_NYSTROEM_COMPONENTS = 300
CLASS_IMBALANCE_TOLERANCE = 0.1  # class_weight="balanced" si |frac_pos-0.5| > esto
SEED = 42


def log(msg):
    print(f"[svm_s12] {msg}", flush=True)


# ---------------------------------------------------------------- features

def build_pair_features(embeddings, pairs, missing_policy="skip"):
    """Devuelve (X, kept_mask). X: hadamard+absdiff [+concat cruda si
    dim<=CONCAT_DIM_THRESHOLD]. kept_mask indica que pares de la lista
    de entrada se conservaron (misma logica de descarte que Seccion 9-11:
    nunca se inventa un embedding faltante)."""
    pairs = list(pairs)
    if not pairs:
        return np.empty((0, 0), dtype=np.float32), np.zeros(0, dtype=bool)
    dim = len(next(iter(embeddings.values())))
    use_concat = dim <= CONCAT_DIM_THRESHOLD
    feats, kept = [], []
    for n1, n2 in pairs:
        e1 = embeddings.get(n1)
        e2 = embeddings.get(n2)
        if e1 is None or e2 is None:
            kept.append(False)
            continue
        parts = [e1 * e2, np.abs(e1 - e2)]
        if use_concat:
            parts.extend([e1, e2])
        feats.append(np.concatenate(parts))
        kept.append(True)
    kept = np.array(kept, dtype=bool)
    if not feats:
        return np.empty((0, 0), dtype=np.float32), kept
    return np.vstack(feats).astype(np.float32), kept


def class_weight_setting(y_train, tol=CLASS_IMBALANCE_TOLERANCE):
    frac_pos = float(np.mean(y_train)) if len(y_train) else 0.5
    setting = "balanced" if abs(frac_pos - 0.5) > tol else None
    return setting, frac_pos


# ---------------------------------------------------------------- modelos

class FittedSVM:
    """Envoltorio comun para el SVM lineal y el RBF-aproximado (Nystroem +
    SGDClassifier hinge), de forma que el resto del pipeline (umbral,
    evaluacion, ranking) no necesite distinguir cual gano en validacion."""

    def __init__(self, kind, clf, sampler=None):
        self.kind = kind  # "linear" | "rbf_approx"
        self.clf = clf
        self.sampler = sampler

    def _transform(self, X):
        return self.sampler.transform(X) if self.sampler is not None else X

    def decision_function(self, X):
        return self.clf.decision_function(self._transform(X))

    def predict_proba(self, X):
        """Adaptador NO probabilistico: solo preserva el orden del score
        continuo (necesario para reutilizar el mismo patron de ranking
        de common_eval.ranked_evaluation, que indexa [:, 1]). No debe
        interpretarse como una probabilidad calibrada."""
        s = self.decision_function(X)
        return np.column_stack([-s, s])


def _make_sgd_svm(c, n_train, cw, seed):
    alpha = 1.0 / (max(c, 1e-6) * max(n_train, 1))
    return SGDClassifier(loss="hinge", alpha=alpha, class_weight=cw, random_state=seed,
                          max_iter=1000, tol=1e-3, early_stopping=True,
                          n_iter_no_change=5, validation_fraction=0.1)


def fit_linear_svm(X_train, y_train, X_val, y_val, c_grid=C_GRID, seed=SEED, cw=None):
    best = None
    for c in c_grid:
        t0 = time.time()
        clf = _make_sgd_svm(c, len(y_train), cw, seed)
        clf.fit(X_train, y_train)
        val_score = clf.decision_function(X_val)
        auprc = average_precision_score(y_val, val_score)
        log(f"    linear C={c}: val_auprc={auprc:.4f} n_iter={clf.n_iter_} ({time.time()-t0:.1f}s)")
        if best is None or auprc > best["auprc"]:
            best = {"clf": clf, "c": c, "auprc": auprc}
    model = FittedSVM("linear", best["clf"])
    diag = {"kernel": "linear", "best_c": best["c"], "val_auprc": best["auprc"]}
    return model, diag


def fit_rbf_approx_svm(X_train, y_train, X_val, y_val, c_grid=C_GRID, gamma_grid=GAMMA_GRID,
                        n_components=N_NYSTROEM_COMPONENTS, seed=SEED, cw=None):
    best = None
    for gamma in gamma_grid:
        t0 = time.time()
        sampler = Nystroem(kernel="rbf", gamma=gamma, n_components=n_components, random_state=seed)
        Xt_train = sampler.fit_transform(X_train)
        Xt_val = sampler.transform(X_val)
        log(f"    rbf gamma={gamma}: Nystroem transform listo ({time.time()-t0:.1f}s)")
        for c in c_grid:
            t1 = time.time()
            clf = _make_sgd_svm(c, len(y_train), cw, seed)
            clf.fit(Xt_train, y_train)
            val_score = clf.decision_function(Xt_val)
            auprc = average_precision_score(y_val, val_score)
            log(f"    rbf gamma={gamma} C={c}: val_auprc={auprc:.4f} n_iter={clf.n_iter_} ({time.time()-t1:.1f}s)")
            if best is None or auprc > best["auprc"]:
                best = {"clf": clf, "sampler": sampler, "c": c, "gamma": gamma, "auprc": auprc}
    model = FittedSVM("rbf_approx", best["clf"], sampler=best["sampler"])
    diag = {
        "kernel": "rbf_nystroem_approx", "best_c": best["c"], "best_gamma": str(best["gamma"]),
        "n_components": n_components, "val_auprc": best["auprc"],
    }
    return model, diag


def select_best_model(X_train, y_train, X_val, y_val, cw=None, seed=SEED):
    """Entrena lineal y RBF-aproximado, se queda con el de mejor AUPRC de
    validacion. Devuelve (modelo_ganador, diagnostico_completo)."""
    lin_model, lin_diag = fit_linear_svm(X_train, y_train, X_val, y_val, cw=cw, seed=seed)
    rbf_model, rbf_diag = fit_rbf_approx_svm(X_train, y_train, X_val, y_val, cw=cw, seed=seed)
    if rbf_diag["val_auprc"] > lin_diag["val_auprc"]:
        winner, winner_diag = rbf_model, rbf_diag
    else:
        winner, winner_diag = lin_model, lin_diag
    full_diag = {
        "linear": lin_diag, "rbf_approx": rbf_diag,
        "elegido": winner_diag["kernel"], "class_weight": cw,
    }
    return winner, full_diag


# ---------------------------------------------------------------- umbral / metricas

def select_threshold(y_val, val_scores):
    """Umbral que maximiza F1 en VALIDACION (nunca en test)."""
    prec, rec, thr = precision_recall_curve(y_val, val_scores)
    if len(thr) == 0:
        return 0.0, {"f1_val": float("nan")}
    f1 = np.zeros(len(thr))
    p, r = prec[:-1], rec[:-1]
    nz = (p + r) > 0
    f1[nz] = 2 * p[nz] * r[nz] / (p[nz] + r[nz])
    best_idx = int(np.argmax(f1))
    return float(thr[best_idx]), {"f1_val": float(f1[best_idx])}


def evaluate_test(y_test, test_scores, threshold):
    y_pred = (test_scores >= threshold).astype(int)
    return {
        "auprc": float(average_precision_score(y_test, test_scores)),
        "auroc": float(roc_auc_score(y_test, test_scores)),
        "f1": float(f1_score(y_test, y_pred, zero_division=0)),
        "precision": float(precision_score(y_test, y_pred, zero_division=0)),
        "recall": float(recall_score(y_test, y_pred, zero_division=0)),
        "n_test": int(len(y_test)),
        "n_pos_test": int(int(y_test.sum())),
        "threshold": float(threshold),
    }


def pr_roc_curves(y_test, test_scores):
    prec, rec, _ = precision_recall_curve(y_test, test_scores)
    fpr, tpr, _ = roc_curve(y_test, test_scores)
    return {
        "precision": prec.astype(np.float32), "recall": rec.astype(np.float32),
        "fpr": fpr.astype(np.float32), "tpr": tpr.astype(np.float32),
    }


# ---------------------------------------------------------------- ranking (Hits@k / MRR)

def ranked_evaluation_svm(model, embeddings, test_positive_pairs, node_pool_t2, rng,
                           k=10, n_negatives_per_query=99, forbidden_pairs=None):
    """Version especifica de Seccion 12 de Hits@k/MRR: reconstruye las
    features candidatas con build_pair_features (misma logica de
    concatenacion condicional que en train/val/test) para que la
    dimensionalidad sea EXACTAMENTE la que el modelo espera -- la
    version generica de common_eval.ranked_evaluation asume siempre
    hadamard+absdiff (256-d), lo cual rompe cuando el modelo se entreno
    con la concatenacion cruda anadida (512-d para fuentes de baja
    dimension). n_negatives_per_query=99 (pool "1 positivo vs 99
    negativos") es el mismo valor por defecto ya usado y documentado en
    common_eval.py; pendiente de confirmacion del usuario antes de
    tratarse como resultado final del TFM (misma nota que en Seccion 9)."""
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
        pairs = [(n1, c) for c in candidates]
        X, kept = build_pair_features(embeddings, pairs, missing_policy="skip")
        if X.shape[0] < 2 or not kept[0]:
            continue
        scores = model.decision_function(X)
        order = np.argsort(-scores)
        # candidates[0] es n2_true; kept preserva el orden relativo, asi
        # que si kept[0]=True (garantizado, ver chequeo arriba) la fila 0
        # de X sigue siendo el positivo verdadero.
        rank = int(np.where(order == 0)[0][0]) + 1
        ranks.append(rank)
    if not ranks:
        return {"hits_at_k": None, "mrr": None, "n_queries": 0, "k": k}
    ranks = np.array(ranks)
    return {
        "hits_at_k": float((ranks <= k).mean()),
        "mrr": float((1.0 / ranks).mean()),
        "n_queries": int(len(ranks)),
        "k": k,
    }
