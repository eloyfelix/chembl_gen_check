# chembl_gen_check

chembl_gen_check is a simple tool that perfoms basic structural sanity checks that can be used as filters or components in a multi-objective optimisation as part of a generative model.

## Installation

```
pip install chembl-gen-check
```

## Usage example

```python
from chembl_gen_check import Checker

checker = Checker()
checker.load_smiles("CC(=O)Oc1ccccc1C(=O)O")

checker.check_scaffold() # Using scaffolds found in ChEMBL
checker.check_ring_systems() # Using ring systems found in ChEMBL
checker.check_lacan() # Profile generated using ChEMBL
checker.check_structural_alerts() # ChEMBL set
```

Code to extract ring systems adapted from: W Patrick Walters. [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils/blob/master/useful_rdkit_utils/ring_systems.py)

Code to calculate LACAN scores adapted from: Dehaen, W. LACAN. https://github.com/dehaenw/lacan/

Using [molbloom](https://github.com/whitead/molbloom) filters for ChEMBL scaffolds and ring systems matching.