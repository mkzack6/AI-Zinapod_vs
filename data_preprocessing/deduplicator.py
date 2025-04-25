import pandas as pd

# Function to resolve SMILES duplicates by keeping the entry with the highest pEC50
def deduplicate_smiles_by_pec50(input_file, output_file, log_file):
    # Load dataset
    df = pd.read_csv(input_file)
    
    # Standardize column names
    df.columns = [col.lower().strip() for col in df.columns]
    required_cols = ["smiles", "ec50", "pec50", "activity"]
    if not all(col in df.columns for col in required_cols):
        raise ValueError("Missing required columns: smiles, ec50, pec50, activity")
    
    # Log initial counts
    initial_count = len(df)
    initial_actives = len(df[df["activity"] == 1])
    initial_decoys = len(df[df["activity"] == 0])
    print(f"Initial dataset: {initial_count} compounds (Actives: {initial_actives}, Decoys: {initial_decoys})")
    
    # Identify duplicates
    duplicate_smiles = df[df["smiles"].duplicated(keep=False)]
    duplicate_count = len(duplicate_smiles)
    unique_duplicate_smiles = duplicate_smiles["smiles"].nunique()
    print(f"Found {duplicate_count} duplicate entries across {unique_duplicate_smiles} unique SMILES")
    
    # Log duplicates for review
    with open(log_file, "w") as f:
        f.write("Duplicate SMILES and their pEC50 values:\n")
        for smiles in duplicate_smiles["smiles"].unique():
            dup_rows = duplicate_smiles[duplicate_smiles["smiles"] == smiles]
            f.write(f"\nSMILES: {smiles}\n")
            for _, row in dup_rows.iterrows():
                f.write(f"  EC50: {row['ec50']}, pEC50: {row['pec50']}, Activity: {row['activity']}\n")
    
    # Sort by SMILES and pEC50 (descending) to keep highest pEC50 first
    df_sorted = df.sort_values(by=["smiles", "pec50"], ascending=[True, False])
    
    # Remove duplicates, keeping the first (highest pEC50)
    df_dedup = df_sorted.drop_duplicates(subset=["smiles"], keep="first")
    
    # Log final counts
    final_count = len(df_dedup)
    final_actives = len(df_dedup[df_dedup["activity"] == 1])
    final_decoys = len(df_dedup[df_dedup["activity"] == 0])
    print(f"After deduplication: {final_count} compounds (Actives: {final_actives}, Decoys: {final_decoys})")
    
    # Save deduplicated dataset
    df_dedup.to_csv(output_file, index=False)
    print(f"Saved deduplicated dataset to {output_file}")
    print(f"Duplicate log saved to {log_file}")

# Main execution
if __name__ == "__main__":
    input_file = "combined_dataset.csv"  # Input CSV with SMILES, EC50, pEC50, Activity
    output_file = "dedup_smiles.csv"     # Output CSV after SMILES deduplication
    log_file = "duplicate_smiles_log.txt"  # Log file for duplicate details
    deduplicate_smiles_by_pec50(input_file, output_file, log_file)