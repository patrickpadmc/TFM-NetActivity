import pandas as pd

df = pd.read_csv("data/processed/analysis/seccion12_svm_link_prediction/seccion12_tabla_completa.tsv", sep="\t",
                  keep_default_na=False)  # para NO tratar el string "None" como NaN

print("--- conteo de rbf_best_gamma (string real, sin re-interpretar 'None' como NaN) ---")
print(df["rbf_best_gamma"].value_counts())

print("\n--- filas con class_weight='balanced' ---")
print(df[df["class_weight"] == "balanced"][["source", "relation", "frac_pos_train", "class_weight"]].to_string(index=False))

print("\n--- filas con modelo_elegido == 'linear' (RBF-aprox NO gano) ---")
print(df[df["modelo_elegido"] == "linear"][["source", "relation", "linear_val_auprc", "rbf_val_auprc"]].to_string(index=False))
