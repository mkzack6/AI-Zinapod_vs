import pandas as pd
import gzip

# Read and extract 1M SMILES
input_file = "CID-SMILES.gz"
output_file = "pubchem_smiles.csv"
max_compounds = 1000000

smiles_list = []
with gzip.open(input_file, "rt") as f:
    for i, line in enumerate(f):
        if i >= max_compounds:
            break
        try:
            cid, smiles = line.strip().split("\t")
            smiles_list.append(smiles)
        except:
            print(f"Skipping invalid line {i+1}: {line.strip()}")
            continue

# Save as CSV
df = pd.DataFrame(smiles_list, columns=["SMILES"])
df.to_csv(output_file, index=False)
print(f"Saved {len(df)} SMILES to {output_file}")