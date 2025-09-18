import pandas as pd
import sys

wk_dir=sys.argv[1]

df = pd.read_csv(f"{wk_dir}/taxonomy/taxonomy.tsv", sep="\t")
# Extract taxonomy, removing species
taxonomy = df["Taxon"]#.str.split(";", expand=True).drop(6, axis=1).values
# taxonomy = pd.Series(taxonomy[:, -1]).str.split("__", expand=True).values[:, -1]
res = pd.Series(taxonomy).value_counts()
res_pct = pd.Series(taxonomy).value_counts(normalize=True)
result = pd.concat([res, res_pct], axis=1)
result.index.name = "classification"
result.to_csv(f"{wk_dir}/taxonomy/taxonomy_table.tsv", sep="\t")
