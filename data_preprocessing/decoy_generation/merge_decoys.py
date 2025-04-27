import pandas as pd
import logging

# Set up logging
logging.basicConfig(
    filename='merge_decoys.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info("Starting script to merge decoys datasets.")

def main():
    # Load the preprocessed decoys
    try:
        preprocessed_df = pd.read_csv('decoys_output_processed_cleaned.csv')
        logging.info(f"Loaded preprocessed decoys with {len(preprocessed_df)} rows.")
        logging.info(f"Columns in preprocessed decoys: {list(preprocessed_df.columns)}")
    except Exception as e:
        logging.error(f"Failed to load decoys_output_processed_cleaned.csv: {e}")
        return

    # Load the DUD-E decoys (no header, so assign column name)
    try:
        dud_e_df = pd.read_csv('dud_e_decoys.csv', header=None, names=['SMILES'])
        logging.info(f"Loaded DUD-E decoys with {len(dud_e_df)} rows.")
        logging.info(f"Assigned column name 'SMILES' to DUD-E decoys data.")
    except Exception as e:
        logging.error(f"Failed to load dud_e_decoys.csv: {e}")
        return

    # Extract Decoy_SMILES from preprocessed decoys
    preprocessed_decoys = preprocessed_df[['Decoy_SMILES']].copy()
    preprocessed_decoys.rename(columns={'Decoy_SMILES': 'SMILES'}, inplace=True)
    preprocessed_decoys['Source'] = 'Preprocessed'

    # Prepare DUD-E decoys
    dud_e_decoys = dud_e_df[['SMILES']].copy()
    dud_e_decoys['Source'] = 'DUD-E'

    # Combine the datasets
    combined_df = pd.concat([preprocessed_decoys, dud_e_decoys], ignore_index=True)
    logging.info(f"Combined dataset has {len(combined_df)} rows before duplicate removal.")

    # Remove duplicates based on SMILES
    initial_rows = len(combined_df)
    combined_df = combined_df.drop_duplicates(subset='SMILES')
    logging.info(f"Removed {initial_rows - len(combined_df)} duplicates. {len(combined_df)} rows remain.")

    # Save the combined dataset
    try:
        combined_df.to_csv('comb_decoys.csv', index=False)
        logging.info("Successfully saved combined decoys to comb_decoys.csv")
    except Exception as e:
        logging.error(f"Failed to save combined decoys: {e}")

if __name__ == "__main__":
    main()