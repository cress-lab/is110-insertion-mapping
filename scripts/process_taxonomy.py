import pandas as pd
import sys

taxonomy_file=sys.argv[1]
output_dir=sys.argv[2]

df = pd.read_csv(taxonomy_file, sep="\t")
# Extract taxonomy, removing species
taxonomy = df["Taxon"]
res = pd.Series(taxonomy).value_counts()
res_pct = pd.Series(taxonomy).value_counts(normalize=True)
result = pd.concat([res, res_pct], axis=1)
result.index.name = "classification"
result.to_csv(f"{output_dir}/taxonomy_table.tsv", sep="\t")
