```markdown
# SMILES Preprocessing for Virtual Screening

This README documents the preprocessing steps applied to the **actives**, **decoys**, and **Zinapod** datasets for virtual screening, focusing on SMILES string standardization and duplicate handling.  
The goal was to ensure chemical validity, consistency, and fairness across all datasets while preserving relevant data for training a robust model.

---

## Overview of Datasets

- **Actives**:  
  A dataset of **1285 compounds** with experimental data (`ec50`, `pec50`, `activity`), used as positive examples in the training set.
  
- **Decoys**:  
  A dataset of **5584 compounds**, used as negative examples in the training set.
  
- **Zinapod**:  
  A screening library of **6774 compounds**, used as the test set for virtual screening.

---

## Preprocessing Steps

The preprocessing pipeline was implemented using **RDKit** and **Python**, with the following steps applied to each dataset:

### SMILES Standardization
- Convert SMILES strings to RDKit molecules.
- Fix nitro groups to correct valence errors (e.g., `[N+](=O)[O-]` → `N(=O)=O`).
- Sanitize molecules, skipping kekulization to avoid altering tautomers unnecessarily.
- Normalize charges, preserving quaternary and aromatic nitrogens.
- Convert back to canonical SMILES using:

  ```python
  Chem.MolToSmiles(canonical=True, isomericSmiles=True, kekuleSmiles=False)
  ```

### Salt Removal and Neutralization
- Remove salts using RDKit’s `SaltRemover`.
- Neutralize charges (except for quaternary/aromatic nitrogens).
- Convert back to canonical SMILES.

### Duplicate Removal
- Identify and remove duplicates based on SMILES strings.
- **Note:** The strategy for duplicate removal differs between datasets (explained below).

---

## Dataset-Specific Preprocessing Results

### Actives Dataset
- **Initial Count**: 1285 rows
- **After Standardization**: 1285 rows (0 invalid SMILES)
- **After Salt Removal/Neutralization**: 1285 rows (no salts removed)
- **Canonical Duplicates**: 503 duplicates identified
- **Duplicate Removal Strategy**:  
  Duplicates were preserved using raw SMILES for duplicate checking.

- **Final Count**: 1285 rows (0% loss)
- **Output File**: `/kaggle/working/actives_processed_preserve_raw.csv`

#### Justification for Keeping Canonical Duplicates
- The **actives** dataset contains experimental metadata (`ec50`, `pec50`, `activity`) per compound.
- Although 503 entries were structurally identical after canonicalization, their metadata may differ due to different experimental conditions.
- Removing these would discard valuable variability critical for robust model training.
- Since **raw SMILES** had 0 duplicates, retaining all 1285 rows preserves all metadata while ensuring chemical validity.

---

### Decoys Dataset
- **Initial Count**: 5584 rows
- **After Standardization**: 5584 rows (0 invalid SMILES)
- **After Salt Removal/Neutralization**: 5584 rows (no salts removed)
- **Canonical Duplicates**: 0 duplicates
- **Duplicate Removal Strategy**:  
  Duplicates were removed based on canonical SMILES (none found).

- **Final Count**: 5584 rows (0% loss)

#### Justification for Removing Canonical Duplicates
- The **decoys** serve as negative examples and lack experimental metadata.
- Structural redundancy would overrepresent certain motifs without adding value.
- Removing duplicates ensures a diverse set ideal for balanced model training.

---

### Zinapod Dataset
- **Initial Count**: 6774 rows
- **After Standardization**: 6765 rows (9 invalid SMILES dropped)
- **After Salt Removal/Neutralization**: 6765 rows (no additional losses)
- **Canonical Duplicates**: 241 duplicates identified
- **Duplicate Removal Strategy**:  
  Duplicates were removed based on canonical SMILES.

- **Final Count**: 6524 rows (3.69% loss)

#### Justification for Removing Canonical Duplicates
- The **Zinapod** dataset is used for virtual screening.
- Structural redundancy increases computational cost without new information.
- Removing duplicates ensures a diverse and efficient screening library.

---

## Why Different Duplicate Strategies?

The preprocessing pipeline aimed to **standardize SMILES strings** across datasets for chemical consistency.  
However, the decision to **keep or remove canonical duplicates** depended on the role and metadata of each dataset:

- **Actives**:  
  Retaining duplicates preserves valuable experimental metadata (e.g., different `ec50` values for the same compound).

- **Decoys and Zinapod**:  
  Lack metadata; thus, canonical duplicates only add redundancy and were removed.

This approach balances chemical standardization with practical requirements for training and screening.

---

## Final Dataset Sizes

| Dataset          | Final Count |
|------------------|-------------|
| Actives          | 1285 rows   |
| Decoys           | 5584 rows   |
| **Training Set** | **6869 rows** |
| Zinapod          | 6524 rows   |

---

## Fairness of the Screening Process

- **Actives**:  
  Retaining 1285 rows ensures dataset diversity, avoiding the 39.14% loss seen with earlier canonical deduplication.
  
- **Decoys**:  
  Structurally unique 5584 compounds create a balanced negative training set.
  
- **Zinapod**:  
  Minimal loss (3.69%) after standardization and deduplication ensures a diverse and efficient screening library.

Overall, each dataset was preprocessed consistently for **chemical validity** while **handling duplicates** in a way that **maximizes the utility** of the data.

---

## How to Use the Processed Data

- **Training Set**:  
  Combine the **actives** (`actives_processed_preserve_raw.csv`) and **decoys** into a single training set (6869 rows total).

- **Screening**:  
  Use the processed **Zinapod** dataset (6524 rows) as the **test set** for virtual screening.

- **Logs**:  
  Refer to log files (e.g., `smiles_preprocessing_preserve_raw.log`) for detailed preprocessing steps and transformations.

---

## Future Considerations

- If structural redundancy in the actives dataset (503 canonical duplicates) becomes a concern (e.g., overrepresentation of motifs), you can rerun preprocessing with canonical deduplication using:

  ```bash
  preprocess_smiles_canonical_dedup.py
  ```

- This would reduce the actives dataset to **782 rows**, but it may also discard valuable experimental metadata.
- **Recommendation**: Retain all 1285 rows unless redundancy significantly impacts model performance.

---
```

---

Would you also like me to help generate a basic `folder structure` diagram (using markdown) for how you'd organize the `/kaggle/working/` directory to accompany this README? 🚀  
It would make the repo even more GitHub-ready!
