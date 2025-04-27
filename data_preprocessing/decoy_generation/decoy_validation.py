from rdkit import Chem
from rdkit.Chem import Descriptors

active_mol = Chem.MolFromSmiles('Cc1ccc(NC(=O)CSc2ncnc3c2cnn3-c2ccccc2Cl)nc1')
decoy_mol = Chem.MolFromSmiles('CC1=CC=CC=C1N2C3=C(C=N2)C(=NC=N3)SCC(=O)NCC4=CC=CS4')

active_mw = Descriptors.MolWt(active_mol)
decoy_mw = Descriptors.MolWt(decoy_mol)
active_logp = Descriptors.MolLogP(active_mol)
decoy_logp = Descriptors.MolLogP(decoy_mol)

print(f"Active MW: {active_mw}, Decoy MW: {decoy_mw}, Diff: {abs(active_mw - decoy_mw)}")
print(f"Active LogP: {active_logp}, Decoy LogP: {decoy_logp}, Diff: {abs(active_logp - decoy_logp)}")