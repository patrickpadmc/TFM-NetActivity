import pandas as pd
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)

df = pd.read_csv("data/processed/analysis/seccion12_svm_link_prediction/seccion12_tabla_completa.tsv", sep="\t")

MODELS_ALL = ["spectral", "node2vec", "hin2vec", "gae", "gat",
              "gae_uniprot_fixed", "gae_uniprot_finetuned",
              "gat_uniprot_fixed", "gat_uniprot_finetuned"]
MODELS_WITH_RAW = MODELS_ALL + ["uniprot_raw", "genept_raw"]

for metric in ["test_auroc", "hits_at_10", "mrr"]:
    piv = df.pivot(index="source", columns="relation", values=metric).reindex(MODELS_WITH_RAW)
    print(f"\n--- {metric} ---")
    print(piv.round(3).to_string())

print("\n--- cobertura (n_train, n_dropped_train, n_test, n_dropped_test) ---")
cov = df[["source", "relation", "n_train", "n_dropped_train", "n_val", "n_dropped_val", "n_test", "n_dropped_test"]]
print(cov.to_string(index=False))

print("\n--- hiperparametros elegidos (modelo, C, gamma, umbral) ---")
hp = df[["source", "relation", "modelo_elegido", "linear_best_c", "rbf_best_c", "rbf_best_gamma",
         "threshold", "f1_val_en_umbral", "class_weight"]]
print(hp.to_string(index=False))
