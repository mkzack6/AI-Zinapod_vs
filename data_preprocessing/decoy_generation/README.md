# Decoy Generation Process

This project outlines the process of generating decoys for a set of active compounds, starting from data preparation to decoy validation. The goal is to produce a dataset of decoys that can be used in virtual screening to evaluate the performance of docking algorithms by distinguishing active compounds from structurally similar but inactive ones.

---

## Overview of the Process

The workflow consists of the following steps:

1. Convert CSV to SMI Format
2. Load and Validate SMILES
3. Compute Morgan Fingerprints
4. Generate Decoys
5. Validate Decoys
6. Preprocess SMILES for Decoys

---

## Step-by-Step Process

### 1. Convert CSV to SMI Format

**Input**: A CSV file containing SMILES strings (e.g., a column of SMILES data).  
**Output**: An SMI file (a simple text file with one SMILES string per line).

**Why**:  
- SMI is a lightweight, widely-used format for SMILES data, compatible with RDKit and other cheminformatics tools.
- Standardizes the input format for downstream processing.

**Process**:

```python
import pandas as pd

# Read CSV
df = pd.read_csv('compounds.csv')

# Write SMILES to SMI file
with open('compounds.smi', 'w') as f:
    for smiles in df['SMILES']:
        f.write(f"{smiles}\n")
```

---

### 2. Load and Validate SMILES

**Input**:  
- `actives.smi` (active compounds)  
- `pubchem_smiles.smi` (decoy pool)

**Output**:  
- Filtered DataFrames with validated SMILES and corresponding RDKit Mol objects.

**Why**:  
- SMILES validation ensures chemically meaningful molecules, avoiding errors in later steps.

**Process**:
- Load SMI files into pandas DataFrames.
- Validate SMILES using RDKit’s `Chem.MolFromSmiles`.
- Filter out invalid SMILES and log issues.

**Dependencies**: `pandas`, `rdkit`

**Criteria**:
- SMILES must parse into valid RDKit Mol objects.
- Molecules with isolated hydrogens are excluded.

---

### 3. Compute Morgan Fingerprints

**Input**: Validated SMILES from actives and decoy pool.  
**Output**: Morgan fingerprints (bit vectors).

**Why**:
- Morgan fingerprints capture substructure information and are ideal for molecular similarity comparison.

**Process**:
- Compute fingerprints using RDKit with `radius=2`, `nBits=2048`.
- Convert fingerprints to NumPy arrays for efficiency.

**Dependencies**: `rdkit`, `numpy`

**Criteria**:
- Fingerprints generated only for validated SMILES.
- Molecules failing sanitization are excluded.

---

### 4. Generate Decoys

**Input**: Morgan fingerprints of actives and decoy pool.  
**Output**: `decoys_output.csv` (active-decoy pairs and Tanimoto similarities).

**Why**:
- Tanimoto similarity is a standard for fingerprint-based molecular similarity.
- A 1:3 active-to-decoy ratio ensures dataset quality and statistical robustness.

**Process**:
- Compute Tanimoto similarity using `sklearn.metrics.pairwise.pairwise_distances` with Jaccard metric.
- Select decoys based on similarity and physicochemical property thresholds.
- Save results to `decoys_output.csv`.

**Dependencies**: `sklearn`, `numpy`, `tqdm`

**Selection Criteria**:
- **Tanimoto Similarity**: High similarity but <1.
- **Physicochemical Properties**:
  - Molecular Weight (MW) difference ≤ 25 Da
  - LogP difference ≤ 1.0
- **Ratio**: 1 active → 3 decoys

---

### 5. Validate Decoys

**Input**: `decoys_output.csv`  
**Output**: Validation logs (`decoy_generation.log`)

**Why**:
- Ensures decoys are high-quality and meet criteria.

**Process**:
- Verify decoy ratio.
- Check Tanimoto similarity range (0.15–0.40 recommended).
- Confirm physicochemical property matches.
- Log statistics and issues.

**Dependencies**: `pandas`, `rdkit`

**Validation Criteria**:
- Each active should have 3 decoys (or fewer if insufficient candidates).
- MW difference ≤ 25 Da.
- LogP difference ≤ 1.0.
- Tanimoto similarities within reasonable range.

---

### 6. SMILES Preprocessing for Decoys Dataset

**Input**: `decoys_output.csv`  
**Output**: `decoys_output_processed_cleaned.csv`

**Why**:
- Further SMILES standardization improves chemical consistency, addressing issues related to charge states, salts, and structural errors.

**Process**:
- Applied `preprocess_smiles.py` to standardize SMILES, remove salts, neutralize charges (except for quaternary and aromatic nitrogens), and eliminate duplicates.
- Corrected common issues:
  - Transformed nitro groups `[N+](=O)[O-]` → `N(=O)=O`.
  - Preserved charges on aromatic nitrogens (e.g., `[n+]`).
  - Skipped kekulization during RDKit sanitization to accommodate charged aromatic systems.

**Dependencies**: `pandas`, `rdkit`

**Results**:

| File | Number of Compounds |
|:-----|:--------------------|
| `decoys_output.csv` | 3855 |
| `decoys_output_processed_cleaned.csv` | 2346 |

- **decrease in Compounds**: The 1:3 ratio has been disturbed.
- **Cleaner Data**: Fewer structural errors and better molecular validity.

**Implications**:
- **Improved Reliability**: The cleaner decoy dataset supports more accurate model training and virtual screening.
- **Decoy-to-Active Ratio**: Adjusted to approximately **1:2.16**, which may affect model balance and needs to be considered during evaluation.
- **Residual Issues**: Some unsanitized SMILES may remain due to the relaxed sanitization settings, though they are minimal.

---

## Dependencies

Install all dependencies using:

```bash
pip install pandas rdkit numpy scikit-learn tqdm
```

> **Note**:  
> Developed and tested with:
> - NumPy: 2.0.2
> - scikit-learn: 1.6.1

---

## Usage

1. **Prepare input data**:
   - Active compounds: `actives.smi`
   - Decoy pool: `pubchem_smiles.smi` (subset of PubChem compounds)

2. **Run the script**:

```bash
python decoy_generation.py
```

3. **Post-process SMILES**:

```bash
python preprocess_smiles.py
```

4. **Check the output**:
   - `decoys_output.csv` — active-decoy pairs with similarity scores
   - `decoys_output_processed_cleaned.csv` — cleaned and standardized SMILES
   - `decoy_generation.log` — logs of validation and issues

---

## Example Output

**Sample `decoys_output.csv`**:

| Active_SMILES | Decoy_SMILES | Tanimoto_Similarity |
| :--- | :--- | :--- |
| Brc1ccc(C(CC2CCCC2)c2cc3cccnc3[nH]2)cc1 | COC1=C(C=C(C=C1)[C@H](CC2=CC=NC=C2)C3=CC=CC=C3)OC4CCCC4 | 0.2535 |

**Validation example**:
- Active MW: 369.31 Da, Decoy MW: 373.50 Da → Difference: 4.19 Da ✅
- Active LogP: 6.04, Decoy LogP: 5.79 → Difference: 0.25 ✅

---

## Limitations and Future Improvements

- **Deprecation Warning**: RDKit’s Morgan fingerprint method shows deprecation warnings; future updates should move to `rdFingerprintGenerator.GetMorganGenerator`.
- **Similarity Range**: Selection could be refined by enforcing a tighter Tanimoto range (e.g., 0.2–0.3).
- **Validation**: More detailed logging of distribution statistics could enhance quality control.
- **Preprocessing Improvements**: Future work could introduce stricter sanitization options or fallback methods for problematic SMILES during preprocessing.

---
