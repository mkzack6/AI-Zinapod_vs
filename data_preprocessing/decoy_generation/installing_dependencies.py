# Suppress output for cleaner logs
import os
os.environ["CONDA_LOG_LEVEL"] = "ERROR"

# Install condacolab to set up Conda in Google Colab
!pip install -q condacolab
import condacolab
condacolab.install()

# Install a newer gcc to update libstdc++ (provides GLIBCXX_3.4.31)
!conda install -c conda-forge gcc=14.1.0 -y

# Install RDKit and other dependencies
!conda install -c conda-forge rdkit=2025.3.1 pandas numpy scikit-learn -y

# Verify libstdc++ version
!strings /usr/local/lib/libstdc++.so | grep GLIBCXX | tail -n 5

# Verify installations
import pkg_resources
required_packages = ['rdkit', 'pandas', 'numpy', 'scikit-learn']
for pkg in required_packages:
    try:
        version = pkg_resources.get_distribution(pkg).version
        print(f"{pkg} version: {version}")
    except pkg_resources.DistributionNotFound:
        print(f"{pkg} is not installed!")

# Import required modules to confirm they work
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, MolStandardize
    from rdkit.Chem.SaltRemover import SaltRemover
    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import MinMaxScaler
    import multiprocessing
    import random
    print("All dependencies imported successfully!")
except ImportError as e:
    print(f"Import failed: {e}")

# Log environment details for debugging
!python --version
!conda list | grep -E 'rdkit|pandas|numpy|scikit-learn|gcc'