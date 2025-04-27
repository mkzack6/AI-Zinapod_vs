import pandas as pd
import logging

# Set up logging
logging.basicConfig(filename='csv_to_smi.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Function to convert CSV to SMI
def csv_to_smi(csv_file, smi_file, smiles_column='smiles'):
    try:
        # Load the CSV file
        df = pd.read_csv(csv_file)
        logging.info(f"Loaded {csv_file} with {len(df)} rows.")

        # Check if the SMILES column exists
        if smiles_column not in df.columns:
            logging.error(f"Column '{smiles_column}' not found in {csv_file}. Available columns: {df.columns}")
            raise ValueError(f"Column '{smiles_column}' not found in {csv_file}.")

        # Extract the SMILES column
        smiles_list = df[smiles_column].dropna().astype(str)
        logging.info(f"Found {len(smiles_list)} SMILES strings in {csv_file}.")

        # Save to SMI file (one SMILES per line)
        with open(smi_file, 'w') as f:
            for smiles in smiles_list:
                f.write(smiles + '\n')
        logging.info(f"Saved {len(smiles_list)} SMILES to {smi_file}.")

    except Exception as e:
        logging.error(f"Failed to convert {csv_file} to {smi_file}: {e}")
        raise

# Convert actives.csv to actives.smi
csv_to_smi('/content/actives.csv', 'actives.smi')

# Convert pubchem.csv to pubchem_smiles.smi
csv_to_smi('/content/pubchem_smiles.csv', 'pubchem_smiles.smi')