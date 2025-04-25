from rdkit import Chem
import pandas as pd

# Load SMILES
df = pd.read_csv("pubchem_smiles.csv")
print(f"Total SMILES: {len(df)}")

# Validate SMILES
invalid_count = 0
for i, smiles in enumerate(df["SMILES"]):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if not mol:
        invalid_count += 1
        print(f"Invalid SMILES at row {i+2}: {smiles}")
    if i % 100000 == 0:
        print(f"Processed {i} SMILES...")

print(f"Invalid SMILES: {invalid_count} ({invalid_count/len(df)*100:.2f}%)")