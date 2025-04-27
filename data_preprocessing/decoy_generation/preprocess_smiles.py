import pandas as pd
import logging
from rdkit import Chem
from rdkit.Chem import SaltRemover

# Set up logging
logging.basicConfig(
    filename='smiles_preprocessing.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info("Starting SMILES preprocessing script.")

def standardize_smiles(smiles):
    """
    Standardize a SMILES string by converting it to a canonical form.
    Returns None if the SMILES is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            logging.warning(f"Invalid SMILES: {smiles}")
            return None
        # Sanitize with permissive settings
        Chem.SanitizeMol(mol, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE))
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logging.warning(f"Failed to standardize SMILES {smiles}: {e}")
        return None

def fix_nitro_groups(mol):
    """
    Convert [N+](=O)[O-] to N(=O)=O to avoid valence errors during sanitization.
    Also handles cases where the nitro group is attached to aromatic rings.
    """
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_CLEANUP))
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N' or atom.GetFormalCharge() != 1:
            continue
            
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
            
        oxygen_atoms = [n for n in neighbors if n.GetSymbol() == 'O']
        if len(oxygen_atoms) != 2:
            continue
            
        o1, o2 = oxygen_atoms
        bond1 = mol.GetBondBetweenAtoms(atom.GetIdx(), o1.GetIdx())
        bond2 = mol.GetBondBetweenAtoms(atom.GetIdx(), o2.GetIdx())
        
        double_bond_o = None
        single_bond_o = None
        for o, bond in [(o1, bond1), (o2, bond2)]:
            if bond.GetBondType() == Chem.BondType.DOUBLE and o.GetFormalCharge() == 0:
                double_bond_o = o
            elif bond.GetBondType() == Chem.BondType.SINGLE and o.GetFormalCharge() == -1 and o.GetDegree() == 1:
                single_bond_o = o
        
        if double_bond_o is None or single_bond_o is None:
            continue
            
        atom.SetFormalCharge(0)
        single_bond_o.SetFormalCharge(0)
        mol.GetBondBetweenAtoms(atom.GetIdx(), single_bond_o.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        logging.info(f"Fixed nitro group in molecule at nitrogen index {atom.GetIdx()}")
    
    return mol

def remove_salts_and_neutralize(smiles):
    """
    Remove salts and neutralize charges in a SMILES string, while preserving necessary charges
    (e.g., quaternary nitrogens and aromatic nitrogens with +1 charge).
    Returns None if processing fails, but attempts to preserve SMILES even with sanitization issues.
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            logging.warning(f"Invalid SMILES for salt removal: {smiles}")
            return None

        # Fix nitro groups
        mol = fix_nitro_groups(mol)

        # Set aromaticity model and sanitize minimally
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        Chem.SanitizeMol(mol, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_CLEANUP))

        # Remove salts
        remover = SaltRemover.SaltRemover()
        mol = remover.StripMol(mol, dontRemoveEverything=True)
        if mol is None or mol.GetNumAtoms() == 0:
            logging.warning(f"No molecule remains after salt removal: {smiles}")
            return None

        # Neutralize charges, but preserve quaternary nitrogens and aromatic nitrogens with +1 charge
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
                # Preserve quaternary nitrogens (4 bonds)
                if atom.GetDegree() == 4:
                    logging.info(f"Preserving charge on quaternary nitrogen at index {atom.GetIdx()} in SMILES {smiles}")
                    continue
                # Preserve aromatic nitrogens with +1 charge (e.g., [n+] in pyridinium-like rings)
                if atom.GetIsAromatic():
                    logging.info(f"Preserving charge on aromatic nitrogen at index {atom.GetIdx()} in SMILES {smiles}")
                    continue
            # Neutralize other charges
            atom.SetFormalCharge(0)

        # Final sanitization, skipping kekulization to avoid errors with charged aromatic nitrogens
        try:
            Chem.SanitizeMol(mol, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES))
        except Exception as e:
            logging.warning(f"Final sanitization failed for SMILES {smiles}, proceeding with unsanitized molecule: {e}")
            return Chem.MolToSmiles(mol, canonical=True)

        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        logging.warning(f"Failed to process SMILES {smiles}: {e}")
        return None

def preprocess_smiles_df(df, smiles_columns):
    """
    Preprocess a DataFrame containing SMILES columns by standardizing, removing salts,
    neutralizing charges, and removing duplicates.
    
    Args:
        df (pd.DataFrame): DataFrame with SMILES columns.
        smiles_columns (list): List of column names containing SMILES strings.
    
    Returns:
        pd.DataFrame: Processed DataFrame.
    """
    logging.info(f"Starting preprocessing with {len(df)} rows.")
    
    for col in smiles_columns:
        logging.info(f"Standardizing SMILES in column: {col}")
        df[f"{col}_standardized"] = df[col].apply(standardize_smiles)
        df = df.dropna(subset=[f"{col}_standardized"])
        logging.info(f"After standardization, {len(df)} rows remain for {col}.")
    
    for col in smiles_columns:
        logging.info(f"Removing salts and neutralizing charges in column: {col}")
        df[f"{col}_processed"] = df[f"{col}_standardized"].apply(remove_salts_and_neutralize)
        df = df.dropna(subset=[f"{col}_processed"])
        logging.info(f"After salt removal and neutralization, {len(df)} rows remain for {col}.")
    
    logging.info("Removing duplicates...")
    df['smiles_tuple'] = df[[f"{col}_processed" for col in smiles_columns]].apply(tuple, axis=1)
    initial_rows = len(df)
    df = df.drop_duplicates(subset='smiles_tuple')
    logging.info(f"Removed {initial_rows - len(df)} duplicates. {len(df)} rows remain.")
    
    for col in smiles_columns:
        df[col] = df[f"{col}_processed"]
    df = df.drop(columns=[f"{col}_standardized" for col in smiles_columns] +
                         [f"{col}_processed" for col in smiles_columns] +
                         ['smiles_tuple'])
    
    return df

def main():
    try:
        df = pd.read_csv('decoys_output.csv')
        logging.info("Successfully loaded decoys_output.csv")
    except Exception as e:
        logging.error(f"Failed to load decoys_output.csv: {e}")
        return
    
    smiles_columns = ['Active_SMILES', 'Decoy_SMILES']
    processed_df = preprocess_smiles_df(df, smiles_columns)
    
    try:
        processed_df.to_csv('decoys_output_processed_cleaned.csv', index=False)
        logging.info("Successfully saved processed data to decoys_output_processed_cleaned.csv")
    except Exception as e:
        logging.error(f"Failed to save processed data: {e}")

if __name__ == "__main__":
    main()