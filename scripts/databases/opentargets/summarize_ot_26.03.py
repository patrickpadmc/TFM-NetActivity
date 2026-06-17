import pyarrow.parquet as pq
import pandas as pd
import glob

files = glob.glob("/home/ppadmoremcc/work/TFM-NetActivity/data/raw/open_targets_26.03/*.parquet")
cols = ["targetId", "diseaseId", "associationScore", "evidenceCount"]

dfs = []
for f in files:
    table = pq.read_table(f, columns=cols)
    dfs.append(table.to_pandas())

df = pd.concat(dfs, ignore_index=True)
print(f"Total filas: {len(df)}")
print(f"Genes únicos: {df['targetId'].nunique()}")
print(f"Enfermedades únicas: {df['diseaseId'].nunique()}")
print("\nDescriptivos de associationScore:")
print(df['associationScore'].describe())
print("\nTop 10 por score:")
print(df.nlargest(10, 'associationScore')[cols].to_string(index=False))
