import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import logging
from sklearn.metrics.pairwise import pairwise_distances
import warnings
from tqdm import tqdm  # For progress bar
warnings.filterwarnings("ignore")  # Suppress RDKit warnings for cleaner output

# Set up logging
logging.basicConfig(filename='decoy_generation.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Load active compounds
try:
    actives_df = pd.read_csv('actives.smi', header=None, names=['SMILES'])
    logging.info(f"Loaded {len(actives_df)} active compounds.")
except Exception as e:
    logging.error(f"Failed to load actives.smi: {e}")
    raise

# Load PubChem SMILES as decoy pool
try:
    decoy_pool_df = pd.read_csv('pubchem_smiles.smi', header=None, names=['SMILES'])
    logging.info(f"Loaded {len(decoy_pool_df)} PubChem SMILES for decoy pool.")
except Exception as e:
    logging.error(f"Failed to load pubchem_smiles.smi: {e}")
    raise

# Function to compute Morgan fingerprint
def get_morgan_fp(smiles, radius=2, nBits=2048):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        return np.array(fp)
    except Exception as e:
        logging.warning(f"Failed to compute Morgan FP for SMILES {smiles}: {e}")
        return None

# Function to check physicochemical properties for decoy selection
def is_valid_decoy(mol, active_mol, mw_diff=25, logp_diff=1.0):
    if mol is None or active_mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    active_mw = Descriptors.MolWt(active_mol)
    logp = Descriptors.MolLogP(mol)
    active_logp = Descriptors.MolLogP(active_mol)
    return (abs(mw - active_mw) <= mw_diff and abs(logp - active_logp) <= logp_diff)

# Compute Morgan fingerprints for actives and decoy pool
actives_df['Mol'] = actives_df['SMILES'].apply(Chem.MolFromSmiles)
actives_df['FP'] = actives_df['SMILES'].apply(get_morgan_fp)
actives_df = actives_df.dropna(subset=['FP'])

decoy_pool_df['Mol'] = decoy_pool_df['SMILES'].apply(Chem.MolFromSmiles)
decoy_pool_df['FP'] = decoy_pool_df['SMILES'].apply(get_morgan_fp)
decoy_pool_df = decoy_pool_df.dropna(subset=['FP'])

logging.info(f"Processed {len(actives_df)} valid active compounds and {len(decoy_pool_df)} valid decoy pool compounds.")

# Select decoys for each active compound
decoys = []
ratio = 3  # 1:3 ratio of actives to decoys
for idx, active_row in tqdm(actives_df.iterrows(), total=len(actives_df), desc="Generating decoys"):
    active_fp = active_row['FP']
    active_mol = active_row['Mol']
    
    # Compute Tanimoto similarity between active and all decoy pool compounds
    decoy_fps = np.array(list(decoy_pool_df['FP']))
    similarities = 1 - pairwise_distances([active_fp], decoy_fps, metric='jaccard')[0]
    
    # Get indices of decoy pool sorted by similarity (descending)
    decoy_indices = np.argsort(similarities)[::-1]
    
    # Select top decoys that meet physicochemical criteria
    selected_decoys = 0
    for decoy_idx in decoy_indices:
        if selected_decoys >= ratio:
            break
        decoy_mol = decoy_pool_df.iloc[decoy_idx]['Mol']
        if is_valid_decoy(decoy_mol, active_mol):
            decoy_smiles = decoy_pool_df.iloc[decoy_idx]['SMILES']
            decoys.append({
                'Active_SMILES': active_row['SMILES'],
                'Decoy_SMILES': decoy_smiles,
                'Tanimoto_Similarity': similarities[decoy_idx]
            })
            selected_decoys += 1

# Save results
decoys_df = pd.DataFrame(decoys)
decoys_df.to_csv('decoys_output.csv', index=False)
logging.info(f"Generated {len(decoys_df)} decoys and saved to decoys_output.csv.")

# Validate the output
if len(decoys_df) == 0:
    logging.warning("No decoys generated. Check SMILES validity or physicochemical criteria.")
else:
    logging.info("Decoy generation completed successfully.")