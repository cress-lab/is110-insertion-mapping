"""
Extract 16S sequences by keyword from the Silva database. 
"""
from Bio import SeqIO

taxonomy_level = -2
to_keep = [
    "Escherichia-Shigella",
    "Klebsiella",
    "Citrobacter",
    # Not sure about these
    "Raoultella",
    "Samsonia",
]
keep_one = False
skip = []

db_path = "/shared/db/silva/latest/fasta/SILVA_138_SSURef_NR99_tax_silva_trunc.replace-whitespace.fasta"
records = SeqIO.parse(db_path, "fasta")

records_desired_taxa = []
for record in records:
    record_taxonomy = record.id.split("_", maxsplit=1)[1].split(";")[taxonomy_level]
    if (record_taxonomy in to_keep) and (record_taxonomy not in skip):
        records_desired_taxa.append(record)
        if keep_one:
            skip.append(record_taxonomy)
            
unique_tag = "_keep_1" if keep_one else ""
print(f"Found {len(records_desired_taxa)} in {', '.join(to_keep)} [level={taxonomy_level}]")

outfile = f"SILVA_138_SSURef_NR99_tax_silva_trunc_subset{len(to_keep)}_level_{taxonomy_level}{unique_tag}.fasta"
with open(f"../sequences/{outfile}", "w") as f:
    SeqIO.write(records_desired_taxa, f, "fasta")