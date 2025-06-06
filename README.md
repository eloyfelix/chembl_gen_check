# chembl_gen_check

chembl_gen_check is a Python library that uses lightweight [MolBloom](https://github.com/whitead/molbloom) filters for rapid verification of the existence of scaffolds, generic scaffolds or ring systems in ChEMBL structures. The library can also indicate whether a compound has uncommon bonds according to the [LACAN](https://github.com/dehaenw/lacan) algorithm, or that a compound triggers a structural alert. Taken together, these checks provide rapid assessment of the reasonableness of ring systems and scaffolds, as well as ensuring that atom and bond environments have precedent.

## Installation

```
pip install chembl-gen-check
```

## Usage example

```python
from chembl_gen_check import Checker

checker = Checker()

smiles = "CCN(CC)C(=O)C[C@H]1C[C@@H]1c1ccccc1"
checker.load_smiles(smiles)

# Murcko scaffold found in ChEMBL (True/False)
checker.check_scaffold()

# Generic Murcko scaffold found in ChEMBL (True/False)
checker.check_skeleton()

# All molecule ring systems found in ChEMBL (True/False)
checker.check_ring_systems()

# Number of structural alerts using the ChEMBL set (integer)
checker.check_structural_alerts()

# LACAN score > 0.5 (True/False)
checker.check_lacan() > 0.5
```

Code to extract ring systems adapted from: W Patrick Walters. [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py)

Code to calculate LACAN scores adapted from: Dehaen, W. LACAN. https://github.com/dehaenw/lacan/
